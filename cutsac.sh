#!/bin/bash
#cut into one hour/ten minute,which means dt is 3600s/600s,pay attention to begin and end
#input parameters:dt(span cutted),delta,(num,which could be canculated by dt/delta),NPTS
while read p1 p2 p3
do
  dt=$p1
  delta=$p2
  NPTS=$p3
done < pmincut
dn=`echo "$dt $delta"|awk '{printf("%i",$1/$2)}'`
   for file in *sac;
   do
      bt=0
      btime=$bt
      et=`echo "$dt $delta"|awk '{printf("%f",$1-$2)}'`
      etime=$et
      et=$dt
      num=$dn
      part=0
#**********************************************************here may be changed if you run this program first time
      station=${file%.*sac}

      for((;num<=NPTS;))
      do
        sac<<EOF
        cut b $bt n $dn
        read $file
        rtrend
        rmean
        hp c 0.1
        lp c 45
        ch b $btime
        ch e $etime
        write ${station}_${part}.SAC
        quit
EOF
      ###########change NPTS in sac to requirement number/fin#########################
      part=$(($part+1))
      bt=`echo "$bt $dt"|awk '{printf("%i",$1+0.25*$2)}'`
#      et=`echo "$bt $dt"|awk '{printf("%.3f",$1+$2)}'`
#      echo $bt $et
      num=$(($num+$dn))
      done
   done
      part=$(($part-1))
      if [ -d cutpart ];then
         rm cutpart
      fi
      echo $part>cutpart

