./LimeScan -f 80M:600M -C 0 -A "LNAW" -w 35M -r 16M -OSR 8 -b 512 -g 48 -n 64 -O TESTfmDabTv -T 1
./LimeScan -f 600M:1000M -C 0 -A "LNAW" -w 35M -r 16M -OSR 8 -b 512 -g 48 -n 64 -O TESTgsmLTE -T 1
./LimeScan -f 1000M:1600M -C 0 -A "LNAW" -w 35M -r 16M -OSR 8 -b 512 -g 48 -n 64 -O TESTsat -T 1
./LimeScan -f 1700M:2400M -C 0 -A "LNAH" -w 35M -r 16M -OSR 8 -b 512 -g 48 -n 64 -O TESTumts -T 1
./LimeScan -f 2400M:2700M -C 0 -A "LNAH" -w 35M -r 16M -OSR 8 -b 512 -g 48 -n 64 -O TESTwifiLTE -T 1
./LimeScan -f 80M:2700M -C 0 -A "LNAW" -w 35M -r 16M -OSR 8 -b 512 -g 48 -n 64  -O TESTall -T 1
#./LimeScan -f 800M:833M -C 0 -A "LNAW" -w 35M -r 16M -OSR 8 -b 256 -g 48.0 -n 16 | tail -5 > log.txt

