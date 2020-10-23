{if ($1 ~/\<START_DATE/) 	{sub($2, START_DATE)}}
{if ($1 ~/\<START_TIME/)	{sub($2, START_TIME)}}
{if ($1 ~/REPORT_START_DATE/) 	{sub($2, START_DATE)}}
{if ($1 ~/REPORT_START_TIME/) 	{sub($2, START_TIME)}}
{if ($1 ~/\<END_DATE/) 		{sub($2, END_DATE)}}
{if ($1 ~/\<END_TIME/) 		{sub($2, END_TIME)}}
{if ($1 ~/\<WET_STEP/) 		{sub($2, WET_STEP)}}
{if ($1 ~/\<REPORT_STEP/) 	{sub($2, REPORT_STEP)}}
{if ($1 ~/\<RAINDAT/) 		{sub($3, substr($3,1,1) RAIN_PATH substr($3, length($3)-5, length($3)))}}
#{if ($1 ~/\<ROUTING_STEP/) 	{sub($2, ROUTING_STEP)}}
{print}

