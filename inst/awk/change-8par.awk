/\[SUBCATCHMENTS\]/ { partof=1 }
/\[SUBAREAS\]/ { partof=0 }
/\[SUBAREAS\]/ { partof4=1 }
/\[INFILTRATION\]/ { partof4=0 }
/\[CONDUITS\]/ { partof6=1 }
/\[XSECTIONS\]/ { partof6=0 }

{ if (partof==1 && $1  ~ /^[0-9]/) { sub ($5, $5*(imperviousness), $5);sub ($6, $6*(width), $6);sub ($7, $7*(slope), $7)}}
{ if (partof4==1 && $1 ~ /^[0-9]/) { sub ($4, $4*simperv, $4);sub ($2, $2*nimperv, $2);sub ($6, $6*(pctzero), $6);sub($5, $5*sperv, $5)}}
{ if (partof6==1 && $1 ~ /^[A-Z0-9]/) { sub ($5, $5*manning, $5)}} # manning = conduits' roughness

{print}
