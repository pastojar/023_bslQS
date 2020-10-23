/<<< Node 288026 >>>/ { partof=1 }
/Link Results/ { partof=0 }
{ if ( partof==1 && $1 ~ /^[JFMASONE]/){print $1 ";" $2 ";" $3}}
{ if ( partof==1 && $1 ~ /DEC/){print $1 ";" $2 ";" $3}}
