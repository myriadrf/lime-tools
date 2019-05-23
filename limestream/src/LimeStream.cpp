/*
Copyright 2019 Evariste COURJAUD - F5OEO
LimeStream is a simple IQ receive or transmit command line

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include <lime/LimeSuite.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(__GNUC__)
#include <unistd.h>
#include <sys/time.h>
#endif
#include <signal.h>
#include <time.h>

// Global variable used by the signal handler and capture/encoding loop
static bool want_quit = false;

float frequency = 434e6;
double bandwidth_calibrating = 2.5e6;
double sample_rate = 2e6;
float gain = 0.8;
unsigned int buffer_size = 1024 * 10;
double postpone_emitting_sec = 0.5;
unsigned int channel = 0;

int mode = -1; // Mode Rx=0,Tx=1,RxTx=2
int buffer_latency = 500;

char *input_filename = nullptr;
char *output_filename = nullptr;

char *calibrationfilename = nullptr;
bool usecalibfile = false;
bool generatecalibfile = false;

FILE *fd_rx = stdout;
FILE *fd_tx = stdin;

FILE *fd_cal = NULL;

int txcount = 0;
int rxcount = 0;

int txfileburst = 1024;
int rxfileburst = 1024;

short *txbuffer[2] = {nullptr, nullptr};
short *rxbuffer[2] = {nullptr, nullptr};

bool verbose = false;

// LimeSDR specific variables
int devicenb = 0; // Index of device to open
lms_device_t *device = NULL;

static void usage()
{

    printf("Usage: \n");
    printf("  -f <Frequency> : Frequency in Hz\n"
           "  -s <Sample rate> : Sample Rate in Hz \n"
           "  -g <Gain (normalized)> : Normalize gain (0.0 to 1.0)\n"
           "  -r <filename> : Receive stream into IQ s16 file (use '-' for stdout)\n"
           "  -t <filename> : Transmit stream from IQ s16 file (use '-' for stdin).\n"
           "  -b <Buffer latency> maximum Buffer latency in milliseconds\n"
           "  -d <Device index> set the Device index if multiple devices\n"
           "  -p <calibration file> : use a Profile calibration file \n"
           "  -c <calibration file> : generate a profile Calibration file \n"
           "  -v : Verbose , print info and statistics \n");
}

static bool validate(lms_device_t *device)
{
    // Get general info about device
    const lms_dev_info_t *device_info;
    device_info = LMS_GetDeviceInfo(device);
    //Get temperature
    float_type Temp;
    LMS_GetChipTemperature(device, 0, &Temp);

    if (device_info != NULL)
    {
    }

    // Get samplerates range
    lms_range_t sr_range[2];
    LMS_GetSampleRateRange(device, LMS_CH_RX, &sr_range[LMS_CH_RX]);
    LMS_GetSampleRateRange(device, LMS_CH_TX, &sr_range[LMS_CH_TX]);

    // Get frequency range
    lms_range_t lo_range[2];
    LMS_GetLOFrequencyRange(device, LMS_CH_RX, &lo_range[LMS_CH_RX]);
    LMS_GetLOFrequencyRange(device, LMS_CH_TX, &lo_range[LMS_CH_TX]);

    if (verbose)
    {
        fprintf(stderr, "%s Library %s Firmware %s Gateware %s ", device_info->deviceName, LMS_GetLibraryVersion(), device_info->firmwareVersion, device_info->gatewareVersion);
        fprintf(stderr, "Temperature %.2f°C\n", Temp);
        fprintf(stderr, "Rx frequency range %.0f-%.0f Mhz\n", lo_range[LMS_CH_RX].min / 1e6, lo_range[LMS_CH_RX].max / 1e6);
        fprintf(stderr, "Tx frequency range %.0f-%.0f Mhz\n", lo_range[LMS_CH_TX].min / 1e6, lo_range[LMS_CH_TX].max / 1e6);
        fprintf(stderr, "Rx samplerate range %.1f-%.1f MSamples/s\n", sr_range[LMS_CH_RX].min / 1e6, sr_range[LMS_CH_RX].max / 1e6);
        fprintf(stderr, "Tx samplerate range %.1f-%.1f MSamples/s\n", sr_range[LMS_CH_TX].min / 1e6, sr_range[LMS_CH_TX].max / 1e6);
    }

    bool checklimit = true;
    int DirectionToCheck = LMS_CH_RX;

    if ((mode == 0) || (mode == 2))
        DirectionToCheck = LMS_CH_RX;
    if ((mode == 1))
        DirectionToCheck = LMS_CH_TX;

    if (sample_rate < sr_range[DirectionToCheck].min)
    {
        fprintf(stderr, "Samplerate too low %f<%f\n",sample_rate,sr_range[DirectionToCheck].min);
        checklimit = false;
    }
    if (sample_rate > sr_range[DirectionToCheck].max)
    {
        fprintf(stderr, "Samplerate too high %f>%f\n",sample_rate,sr_range[DirectionToCheck].max);
        checklimit = false;
    }
    if (frequency < lo_range[DirectionToCheck].min)
    {
        fprintf(stderr, "Frequency too low %f<%f\n",frequency,lo_range[DirectionToCheck].min);
        checklimit = false;
    }
    if (frequency > lo_range[DirectionToCheck].max)
    {
        fprintf(stderr, "Frequency too high %f>%f\n",frequency,lo_range[DirectionToCheck].max);
        checklimit = false;
    }

    return checklimit;
}

#ifdef _MSC_VER
BOOL WINAPI
signal_handler(int signum)
{
    if (CTRL_C_EVENT == signum)
    {
        fprintf(stdout, "Caught signal %d\n", signum);
        want_quit = true;
        return TRUE;
    }
    return FALSE;
}
#else
void signal_handler(int signum)
{
    fprintf(stdout, "Caught signal %d\n", signum);
    want_quit = true;
}
#endif

int main(int argc, char **argv)
{
    int opt;
    bool result = false;

    while ((opt = getopt(argc, argv, "f:s:g:r:t:b:d:p:c:vh?")) != EOF)
    {
        result = true;
        switch (opt)
        {
        case 'f':
            frequency = atof(optarg);
            break;
        case 's':
            sample_rate = atof(optarg);
            break;
        case 'g':
            gain = atof(optarg);
            break;
        case 'r':
            output_filename = optarg;
            mode=0; //rx
            break;
        case 't':
            input_filename = optarg;
            mode=1; //tx
            break;
        case 'b':
            buffer_latency = atoi(optarg);
            break;
        case 'd':
            devicenb = atoi(optarg);
            break;
        case 'p':
            calibrationfilename = optarg;
            usecalibfile = true;
            break;
        case 'c':
            calibrationfilename = optarg;
            generatecalibfile = true;
            break;
        case 'v':
            verbose = true;
            break;
        case 'h':
        case '?':
            usage();
            return EXIT_SUCCESS;

        default:
            fprintf(stderr, "unknown argument '-%c %s'\n", opt, optarg);
            usage();
            return EXIT_FAILURE;
        }

        if (!result)
        {
            fprintf(stderr, "argument error: '-%c %s'\n", opt, optarg);
            usage();
            return EXIT_FAILURE;
        }
    }

    if (!result)
    {
        usage();
        return EXIT_FAILURE;
    }

    if(mode==-1)
    {
        fprintf(stderr,"Specify at least -t for transmit or -r to receive \n");
        usage();
        return EXIT_FAILURE;
    }

    bandwidth_calibrating = sample_rate;

    //Calculate buffersize regarding samplerate and latency
    // Fpga packet is is 1020 * IQ16 samples (USB burst)
    //double rate = sample_rate/1e6;
    //streamSize = (mTxStreams[0].used||mRxStreams[0].used) + (mTxStreams[1].used||mRxStreams[1].used);
    //rate = (rate + 5) * config.performanceLatency * streamSize;
    // rate is multiple of 2
    // For 4MS one channel , (4+5)=9 mulitple of 2^n -> 8
    // Fifo is performing 8*1020 samples at once

    double rate = sample_rate / 1e6;
    rate += 5;
    int BatchSize = 0;
    for (int batch = 1; batch < rate; batch <<= 1)
    {
        BatchSize = batch;
    }

    float SampleDelaySize = sample_rate * buffer_latency * 1e-3;
    txfileburst = BatchSize * 1020;
    rxfileburst = BatchSize * 1020;
    if (SampleDelaySize < BatchSize * 1020)
    {
        SampleDelaySize = BatchSize * 1020;

        fprintf(stderr, "Latency too short, setting to %f\n", SampleDelaySize);
    }
    buffer_size = (int)SampleDelaySize;
    buffer_latency = SampleDelaySize * 1e3 / sample_rate;
    int burstlatency = ((BatchSize * 1020) * 1e3) / sample_rate;
    fprintf(stderr, "Info:Buffer size=%d samples, Buffer latency=%d ms, Burst latency=%d ms\n", buffer_size, buffer_latency, burstlatency);

    // Let's catch some signals
#ifdef _MSC_VER
    SetConsoleCtrlHandler((PHANDLER_ROUTINE)signal_handler, TRUE);
#else
    signal(SIGINT, &signal_handler);
    signal(SIGILL, &signal_handler);
    signal(SIGFPE, &signal_handler);
    //signal(SIGSEGV, &signal_handler);
    signal(SIGTERM, &signal_handler);
    signal(SIGABRT, &signal_handler);
#endif

    // Open file for RX
    if ((mode == 0) || (mode == 2)) //Rx
    {
        if (strcmp(output_filename, "-") == 0)
        {
            fd_rx = stdout;
        }
        else
        {
            fd_rx = fopen(output_filename, "wb");
        }
        if (fd_rx == NULL)
        {
            fprintf(stderr, "Failed to open file: %s\n", output_filename);
            return EXIT_FAILURE;
        }

        result = setvbuf(fd_rx, NULL, _IOFBF, rxfileburst);
        if (result != 0)
        {
            fprintf(stderr, "setvbuf() failed: %d\n", result);
            return EXIT_FAILURE;
        }

        rxcount++;
    }

    // Open file for TX
    if ((mode == 1) || (mode == 2)) //Tx
    {
        if (strcmp(input_filename, "-") == 0)
        {
            fd_tx = stdin;
        }
        else
        {
            fd_tx = fopen(input_filename, "rb");
        }
        if (fd_tx == NULL)
        {
            fprintf(stderr, "Failed to open file: %s\n", input_filename);
            return EXIT_FAILURE;
        }

        result = setvbuf(fd_tx, NULL, _IOFBF, txfileburst);
        if (result != 0)
        {
            fprintf(stderr, "setvbuf() failed: %d\n", result);
            return EXIT_FAILURE;
        }
        txcount++;
    }

    fprintf(stderr, "Alloc Buffer Tx %d, Rx %d \n", txcount, rxcount);

    // Allocation of intermediate I16 buffers
    for (int i = 0; i < txcount; i++)
    {
        txbuffer[i] = (short *)malloc(txfileburst * sizeof(short) * 2);
    }
    for (int i = 0; i < rxcount; i++)
    {
        rxbuffer[i] = (short *)malloc(rxfileburst * sizeof(short) * 2);
    }

// ******* Starting working with Lime *************

// Detect if a Lime device is there
#define MAX_DEVICE 8
    lms_info_str_t device_list[MAX_DEVICE];
    int devicecount = 0;
    if ((devicecount = LMS_GetDeviceList(device_list)) < 0)
    {
        fprintf(stderr, "LMS_GetDeviceList() : %s\n", LMS_GetLastErrorMessage());
        return EXIT_FAILURE;
    }

    // Choose the the device n
    if (devicenb < devicecount)
    {
        LMS_Open(&device, device_list[devicenb], NULL);
    }
    else
    {
        fprintf(stderr, "Device : %d not present(should be <%d)\n", devicenb, devicecount - 1);
        return EXIT_FAILURE;
    }

    //Init device
    if (LMS_Init(device) < 0)
    {
        fprintf(stderr, "LMS_Init() : %s\n", LMS_GetLastErrorMessage());
        return EXIT_FAILURE;
    }

    // Validate parameters
    if (!validate(device))
        return EXIT_FAILURE;

    // Enable tx channels
    for (int i = 0; i < txcount; i++)
    {
        if (LMS_EnableChannel(device, LMS_CH_TX, i, true) < 0)
        {
            fprintf(stderr, "LMS_EnableChannelTx() : %s\n", LMS_GetLastErrorMessage());
            return EXIT_FAILURE;
        }
        // Avoid having a spike at start
        LMS_SetNormalizedGain(device, LMS_CH_TX, i, 0);
    }

    // Enable rx channels
    for (int i = 0; i < rxcount; i++)
    {
        if (LMS_EnableChannel(device, LMS_CH_RX, i, true) < 0)
        {
            fprintf(stderr, "LMS_EnableChannelRx() : %s\n", LMS_GetLastErrorMessage());
            return EXIT_FAILURE;
        }
        LMS_SetNormalizedGain(device, LMS_CH_RX, i, 0);
    }

    // Set common samplerate rx/tx and let limesuite decide upsample/downsample
    if (LMS_SetSampleRate(device, sample_rate, 0) < 0)
    {
        fprintf(stderr, "LMS_SetSampleRate() : %s\n", LMS_GetLastErrorMessage());
        return EXIT_FAILURE;
    }

    // Set rx frequency : Fixme! Need several frequencies parameters
    for (int i = 0; i < rxcount; i++)
    {
        if (LMS_SetLOFrequency(device, LMS_CH_RX, i, frequency) < 0)
        {
            fprintf(stderr, "LMS_SetLOFrequency() : %s\n", LMS_GetLastErrorMessage());
            return EXIT_FAILURE;
        }
    }

    // Set tx frequency : Fixme! Need several frequencies parameters
    for (int i = 0; i < txcount; i++)
    {
        if (LMS_SetLOFrequency(device, LMS_CH_TX, i, frequency) < 0)
        {
            fprintf(stderr, "LMS_SetLOFrequency() : %s\n", LMS_GetLastErrorMessage());
            return EXIT_FAILURE;
        }
    }

    // Set rx gain : Fixme! Need several gain parameters for each
    for (int i = 0; i < rxcount; i++)
    {
        if (LMS_SetNormalizedGain(device, LMS_CH_RX, i, gain) < 0)
        {
            fprintf(stderr, "LMS_SetNormalizedGain() : %s\n", LMS_GetLastErrorMessage());
            return EXIT_FAILURE;
        }
    }

    // Set tx gain : Fixme! Need several gain parameters for each
    for (int i = 0; i < txcount; i++)
    {
        if (LMS_SetNormalizedGain(device, LMS_CH_TX, i, gain) < 0)
        {
            fprintf(stderr, "LMS_SetNormalizedGain() : %s\n", LMS_GetLastErrorMessage());
            return EXIT_FAILURE;
        }
    }

    // Calibration process : warning, if a channel is not enable, calibration is not working
    if (!usecalibfile)
    {

        // Perform Rx Calibration using LimeSuite calibration process
        for (int i = 0; i < rxcount; i++)
        {

            if (LMS_Calibrate(device, LMS_CH_RX, i, bandwidth_calibrating, 0) < 0)
            {
                fprintf(stderr, "LMS_Calibrate() : %s\n", LMS_GetLastErrorMessage());
                return EXIT_FAILURE;
            }
        }

        // Perform Tx Calibration using LimeSuite calibration process : Warning , Spikes output
        for (int i = 0; i < txcount; i++)
        {

            if (LMS_Calibrate(device, LMS_CH_TX, i, bandwidth_calibrating, 0) < 0)
            {
                fprintf(stderr, "LMS_Calibrate() : %s\n", LMS_GetLastErrorMessage());
                return EXIT_FAILURE;
            }
        }

        if (generatecalibfile)
        {
            // Save a calibration file
            if (LMS_SaveConfig(device, calibrationfilename) < 0)
            {
                fprintf(stderr, "LMS_SaveConfig() : %s\n", LMS_GetLastErrorMessage());
                return EXIT_FAILURE;
            }
        }
    }
    else
    {
        // Load from a calibration file : warning, all parameters overidden : FixMe should be only IQ offset parameter
        if (LMS_LoadConfig(device, calibrationfilename) < 0)
        {
            fprintf(stderr, "LMS_LoadConfig() : %s\n", LMS_GetLastErrorMessage());
            return EXIT_FAILURE;
        }
    }

    // Configure Rx stream
    lms_stream_t rx_stream[2];
    for (int i = 0; i < rxcount; i++)
    {

        rx_stream[i].channel = i;
        rx_stream[i].fifoSize = buffer_size;
        rx_stream[i].throughputVsLatency = 1;
        rx_stream[i].isTx = LMS_CH_RX;
        rx_stream[i].dataFmt = lms_stream_t::LMS_FMT_I16;

        if (LMS_SetupStream(device, &rx_stream[i]) < 0)
        {
            fprintf(stderr, "LMS_SetupStream() : %s\n", LMS_GetLastErrorMessage());
            return EXIT_FAILURE;
        }
        if (LMS_EnableChannel(device, LMS_CH_RX, i, true) < 0)
        {
            fprintf(stderr, "LMS_EnableChannelRx() : %s\n", LMS_GetLastErrorMessage());
            return EXIT_FAILURE;
        }
        LMS_StartStream(&rx_stream[i]);
    }

    // Configure Tx stream
    lms_stream_t tx_stream[2];
    for (int i = 0; i < txcount; i++)
    {

        tx_stream[i].channel = i;
        tx_stream[i].fifoSize = buffer_size;
        tx_stream[i].throughputVsLatency = 1;
        tx_stream[i].isTx = LMS_CH_TX;
        tx_stream[i].dataFmt = lms_stream_t::LMS_FMT_I16;
        if (LMS_SetupStream(device, &tx_stream[i]) < 0)
        {
            fprintf(stderr, "LMS_SetupStream() : %s\n", LMS_GetLastErrorMessage());
            return EXIT_FAILURE;
        }
        if (LMS_EnableChannel(device, LMS_CH_TX, i, true) < 0)
        {
            fprintf(stderr, "LMS_EnableChannelTx() : %s\n", LMS_GetLastErrorMessage());
            return EXIT_FAILURE;
        }
        LMS_StartStream(&tx_stream[i]);
    }

    // Main loop for streaming

    /*bool rx_isapipe=false;
    if(rxcount>0)
         rx_isapipe= (fseek(fd_rx, 0, SEEK_CUR) < 0); //Dirty trick to see if it is a pipe or not
    */
    bool tx_isapipe = false;
    if (txcount > 1)
        tx_isapipe = (fseek(fd_tx, 0, SEEK_CUR) < 0); //Dirty trick to see if it is a pipe or not

    int samplecount = 0;
    lms_stream_status_t status;

    while (!want_quit)
    {
        //Rx
        for (int i = 0; i < rxcount; i++)
        {
            int nb_samples_received = LMS_RecvStream(&rx_stream[i], rxbuffer[i], rxfileburst, nullptr, buffer_latency);
            samplecount += nb_samples_received;
            if (nb_samples_received != rxfileburst)
            {
                fprintf(stderr, "Receiving timeout, only %d/%d bytes received\n", nb_samples_received, rxfileburst);
            }
            int nb_samples_written = fwrite(rxbuffer[i], sizeof(rxbuffer[i]), rxfileburst, fd_rx);
            if (nb_samples_written != rxfileburst)
            {
                // Pipe or regular file : not enough space, quit immediatelly
                fprintf(stderr, "Write only %d bytes on output\n", nb_samples_written);
                want_quit = true;
                break;
            }

            LMS_GetStreamStatus(&rx_stream[i], &status); //Get stream status
        }

        //Tx
        for (int i = 0; i < txcount; i++)
        {
            int nb_samples_to_send = fread(txbuffer[i], sizeof(txbuffer[i]), txfileburst, fd_tx);
            if (nb_samples_to_send != txfileburst)
            {
                if (!tx_isapipe)
                {
                    //We end with a regular file, we quit but send the last packet
                    want_quit = true;
                }
                else
                {
                    fprintf(stderr, "Read only %d bytes on input\n", nb_samples_to_send);
                }
            }
            if (nb_samples_to_send > 0)
            {
                int nb_samples_sent = LMS_SendStream(&tx_stream[i], txbuffer[i], nb_samples_to_send, nullptr, buffer_latency);
                if (nb_samples_sent != txfileburst)
                {
                    fprintf(stderr, "Sending timeout, only %d/%d bytes sent\n", nb_samples_sent, txfileburst);
                }
                samplecount += nb_samples_sent;
            }

            LMS_GetStreamStatus(&tx_stream[i], &status); //Get stream status
        }
        if (verbose)
        {
            if (samplecount >= (int)sample_rate)
            {
                fprintf(stderr, "Status : \tSampleRate %.0f \tfifo %d%% \tunderflow %d \toverflows %d \tdrop %d\n", status.linkRate / 4.0, (status.fifoFilledCount * 100) / status.fifoSize, status.underrun, status.overrun, status.droppedPackets);
                samplecount = 0;
            }
        }
    }
    //Wait at least on 2* buffer latency to send the latest samples if any
    usleep(buffer_latency * 1e3);

    //Close device
    LMS_Close(device);

    // Free of intermediate I16 buffers

    for (int i = 0; i < txcount; i++)
    {
        free(txbuffer[i]);
    }
    for (int i = 0; i < rxcount; i++)
    {
        free(rxbuffer[i]);
    }
}
