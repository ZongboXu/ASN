#!/bin/bash

#touch time
#echo "`date`">time
NPTS=960000
#hour
dt=60
delta=0.01
#station begin/end station
bs=25
es=25
#station for cross-correlation
SPATH="/home/zongboxu/Documents"
###############
dn=`echo "$dt $delta"|awk '{printf("%i",$1/$2)}'`
###############
cd $SPATH/BSU/EVEN/SPACE/
cp $SPATH/PROGRAM/bash/cutsac.sh  ./cutsac.sh
#cp $SPATH/PROGRAM/bash/timenorm.sh ./timenorm.sh
chmod 771 cutsac.sh
break
###############
#list all station names in arrstation
nstation=1
for file in *.sac;
do
   exist=0
#*********need be changed if you run this program first time
#   station=${file%%_*}
   station=${file%.*}
#**********
   if [ "$nstation" == "1" ];then     
   arrstation[$nstation]=$station
   nstation=$(($nstation+1))
   else
    for((i=1;i<=$(($nstation-1));i++ ))
    do
       if [ "${arrstation[$i]}" == "$station" ];then
          exist=1
          break
       fi
    done
       if [ "$exist" == "1" ];then
          continue
       else

       arrstation[$nstation]=$station
       nstation=$(($nstation+1))
       fi
   fi
done
echo ${arrstation[$bs]},$nstation
#the number of station name 
nstation=$(($nstation-1))
#exit 
##################cut and timenormalization
      if [ -d pmincut ];then
         rm pmincut
      fi
      touch pmincut
      echo $dt" "$delta" "$NPTS>pmincut

sh cutsac.sh
##sh timenorm.sh
###############array of time *******#############################
ntstation=1
#assuming that there are only ONE type of stations' names and collect all time name in arrtstation  
for file in *.SAC;
do
   exist=0
#******************need be changed if you run this program first time
   tstation=${file#*_}
#**********************
   if [ "$ntstation" == "1" ];then     
   arrtstation[$ntstation]=$tstation
   ntstation=$(($ntstation+1))
   else
    for((i=1;i<=$(($ntstation-1));i++ ))
    do
       if [ "${arrtstation[$i]}" == "$tstation" ];then
          exist=1
          break
       fi
    done
       if [ "$exist" == "1" ];then
          continue
       else

       arrtstation[$ntstation]=$tstation
       ntstation=$(($ntstation+1))
       fi
   fi
done
#the number of station name 
ntstation=$(($ntstation-1))
#################cross-correlating###############################
for((;bs<=$es;bs=$(($bs+1))))
do
  for((pf=1;pf<=$nstation;pf=$(($pf+1))))
  do
    for((n=1;n<=$ntstation;))
    do
      sac <<EOF
      read ${arrstation[$bs]}_${arrtstation[n]} ${arrstation[$pf]}_${arrtstation[n]}
      correlate master 2
      write ${arrstation[$bs]}${arrstation[$pf]}_${arrtstation[$n]} tempfile
      quit
EOF
      n=$(($n+1))
    done
    for tfile in ${arrstation[$bs]}${arrstation[$pf]}_*;
    do
      sac <<EOF
      read $tfile
      write tempfile
      quit
EOF
    break
    done
    for file in ${arrstation[$bs]}${arrstation[$pf]}_*;
    do
      sac <<EOF
      read tempfile
      addf $file
      write tempfile
      quit
EOF
    done
      sac <<EOF
      read tempfile
      subf $tfile
      write tempfile
      quit
EOF
    echo ${arrstation[$bs]} AND ${arrstation[$pf]} stackdone
    mv ./tempfile ./${arrstation[$bs]}${arrstation[$pf]}.cor
  done
done
#cd ~/ABSN/bash
#echo "`date`">>time
exit 0
