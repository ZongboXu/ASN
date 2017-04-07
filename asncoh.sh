#!/bin/bash

#touch time
#echo "`date`">time
NPTS=960000
#hour
dt=60
delta=0.01
#station begin/end station
bs=25
#25
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

#sh cutsac.sh
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
#################cross-coherence###############################
     if [ -d snamelist ];then
         rm snamelist
     fi
     touch snamelist
     echo ${arrstation[*]}>snamelist
     if [ -d snamenum ];then
         rm snamenum
     fi
     touch snamenum
     echo $nstation>snamenum
     if [ -d tnamelist ];then
         rm tnamelist
     fi
     touch tnamelist
     echo ${arrtstation[*]}>tnamelist
     if [ -d tnamenum ];then
         rm tnamenum
     fi
     touch tnamenum
     echo $ntstation>tnamenum
     if [ -d loop ];then
         rm loop
     fi
     touch loop
     echo $bs $es>loop

$SPATH/PROGRAM/fortran/Coherence/asncoh
#$SPATH/PROGRAM/fortran/Coherence/SAasncoh
#cd ~/ABSN/bash
#echo "`date`">>time
exit 0
