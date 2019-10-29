#
# Copyright 2019 Lime Microsystems Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#./LTEscan -B 1,3,5,7,8,20,28 -b 512 -A LNAW
./LTEscan -B 20 -n 2 -A LNAL
# follows convention of LimeScan software where possible
# -b NFFT size {128,256,512,1024,512} default 512
# List of LTE bands, separated by commas and no spaces e.g. 1,5,7,20  Not all bands may be available in your area.
# -A antenna port {LNAH,LNAW,LNAL} available channels depend on board type.  LNAW for LimeSDR-Mini, LNAL for LimeNet Micro.
# -C channel number {0,1}, default 0
# -f frequency correction in Hz.  Typically between -800Hz and 800Hz, correcting for TXCO error.
# -v --version
# -h --help
# --info not implemented
# -g channel gain 48 (Do not alter)
# -OSR over sample {1,2,4,8,...} Max value depends on NFFT, 2048 OSR<=4, 1024 OSR<=8, 512 OSR<=16 256 OSR<=32 128 OSR<=32
# -o -O output path/filename stub (no extention, currently not used)
#
# Results in LTEscan.html and ./output/LTEscan.txt
# Detailed diagnostics in ./output/LTElog.txt
