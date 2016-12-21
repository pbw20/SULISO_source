#!/usr/bin/awk
# this awk script parses the .fort.7 file for coordinates
BEGIN { 
	elem[1] = "  1H "
	elem[6] = " 12C "
	elem[7] = " 14N "
	elem[8] = " 16O "
	elem[9] = "   F "
	elem[17] = "  CL "
	elem[35] = "  BR "
	elem[44] = "  RU "
	}
{
	if ( NF == 4 ){
	gsub(/\D/,"E")
	printf( "%s%20.10E%20.10E%20.10E\n",elem[$1],$2,$3,$4)}
	}
{
} 
