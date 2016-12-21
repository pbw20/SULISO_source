#!/bin/sh
# $2 is the internal coordinate file, $1 is the .fort.7 file, $3 is the .vib file
if [[ -f $2 && -f $1 ]]
then
	cat ~/bin/SULISO/vib_header.txt > $3
	awk -f ~/bin/SULISO/process_crd.awk $1 >> $3
	cat ~/bin/awks/symm_info_6_combo.txt >>$3
	numheaderlines=`wc -l ~/bin/SULISO/vib_header.txt | awk '{print $1}'`
	numoutlines=`wc -l $3 | awk '{print $1}'`
	numatoms=`expr $numoutlines - $numheaderlines`

	awk -f ~/bin/SULISO/process_hess.awk $1 >> $3
	dos2unix $3
	~/bin/SULISO/camvib2016.exe < $3 >vibout.txt
else
	echo "Usage: ./run_camvib <fort.7 file> <symm_info file> <output file (vib)>"
fi
exit 0
