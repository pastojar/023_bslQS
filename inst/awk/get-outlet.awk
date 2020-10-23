/<<< Node 288026 >>>/ { partof=1 }
/Link Results/ { partof=0 }
{ if ( partof==1 && $1 ~ /^[JFMASON]/){print $3 ";" $2} }
{ if ( (partof==1) && ($1 ~ /DEC/) )  {print $3 ";" $2} }