#!/usr/bin/awk
# this awk script parses the fort.7 file for the Hessian
BEGIN { skip_count=0 }
#
#{
#if ( NR == 1 ) { 
#	next }
#}

{
if ( NF == 4 ) {
  skip_count++
  next
  } 
}

{
#	if (NR > (2 * skip_count) + 1)
	if (NR > (2 * skip_count)) {
	gsub(/\D/,"E")
	print $0}
#	if (NR > (skip_count)) print $0
#	if (NR > 8) print $0
#	if (NR > 106) print $0
#	if (NR > 26) print $0
}

END { print "STOP" }
