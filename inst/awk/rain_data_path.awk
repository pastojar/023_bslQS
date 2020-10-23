{if ($1 ~/\<RAINDAT/) 		{print(RAIN_PATH "/" substr($3,2, length($3)-5 ) "dat")}}

