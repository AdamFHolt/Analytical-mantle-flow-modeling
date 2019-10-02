#!/bin/bash

# global parameters
smooth=0
filter=1
box_width=30 	# only for smooth = 1
dip_int=250.0

# south to north
lat_min=12.5; lat_max=35; dip_depth=-300.0000; dip_depth_shall=-260.0000; dip_dir=0;
slab2_cont='izu_slab2_dep_02.24.18_contours.in'; out1='Izu';
./extract_dip_north-south_2Contours $lat_min $lat_max $dip_depth $dip_depth_shall $dip_dir $smooth $box_width $dip_int $slab2_cont $out1

lat_min=23.5; lat_max=32.5; dip_depth=-300.0000; dip_depth_shall=-260.0000; dip_dir=0;
slab2_cont='ryu_slab2_dep_02.26.18_contours.in'; out7='Ryu';
./extract_dip_north-south_2Contours $lat_min $lat_max $dip_depth $dip_depth_shall $dip_dir $smooth $box_width  $dip_int $slab2_cont $out7

lat_min=35; lat_max=55; dip_depth=-300.0000; dip_depth_shall=-260.0000; dip_dir=0;
slab2_cont='kur_slab2_dep_02.24.18_contours.in'; out2='Jap';
./extract_dip_north-south_2Contours $lat_min $lat_max $dip_depth $dip_depth_shall $dip_dir $smooth $box_width  $dip_int $slab2_cont $out2

lat_min=-39; lat_max=-19; dip_depth=-300.0000; dip_depth_shall=-260.0000; dip_dir=0;
slab2_cont='ker_slab2_dep_02.24.18_contours.in'; out3='Ton';
./extract_dip_north-south_2Contours $lat_min $lat_max $dip_depth $dip_depth_shall $dip_dir $smooth $box_width  $dip_int $slab2_cont $out3

lat_min=-44; lat_max=7; dip_depth=-300.0000; dip_depth_shall=-260.0000; dip_dir=1;
slab2_cont='sam_slab2_dep_02.23.18_contours.in'; out4='SAm';
./extract_dip_north-south_2Contours $lat_min $lat_max $dip_depth $dip_depth_shall $dip_dir $smooth $box_width  $dip_int $slab2_cont $out4

lat_min=7; lat_max=20; dip_depth=-240.0000; dip_depth_shall=-200.0000; dip_dir=1;
slab2_cont='cam_slab2_dep_02.24.18_contours.in'; out5='CAm';
./extract_dip_north-south_2Contours $lat_min $lat_max $dip_depth $dip_depth_shall $dip_dir $smooth $box_width  $dip_int $slab2_cont $out5

lat_min=39; lat_max=49; dip_depth=-260.0000; dip_depth_shall=-220.0000; dip_dir=1;
slab2_cont='cas_slab2_dep_02.24.18_contours.in'; out6='NAm';
./extract_dip_north-south_2Contours $lat_min $lat_max $dip_depth $dip_depth_shall $dip_dir $smooth $box_width  $dip_int $slab2_cont $out6

lat_min=-60.2; lat_max=-55.1; dip_depth=-260.0000; dip_depth_shall=-220.0000; dip_dir=0;
slab2_cont='sco_slab2_dep_02.23.18_contours.in'; out8='Scot';
./extract_dip_north-south_2Contours $lat_min $lat_max $dip_depth $dip_depth_shall $dip_dir $smooth $box_width  $dip_int $slab2_cont $out8

lat_min=-20.5; lat_max=-11.5; dip_depth=-260.0000; dip_depth_shall=-220.0000; dip_dir=1;
slab2_cont='van_slab2_dep_02.23.18_contours.in'; out9='Van';
./extract_dip_north-south_2Contours $lat_min $lat_max $dip_depth $dip_depth_shall $dip_dir $smooth $box_width  $dip_int $slab2_cont $out9

# west to east
lon_min=90; lon_max=118; dip_depth=-260.0000; dip_depth_shall=-220.0000; dip_dir=0;
slab2_cont='sum_slab2_dep_02.23.18_contours.in'; out10='Java';
./extract_dip_east-west_2Contours $lon_min $lon_max $dip_depth $dip_depth_shall $dip_dir $smooth $box_width  $dip_int $slab2_cont $out10

lon_min=178; lon_max=206.5; dip_depth=-200.0000; dip_depth_shall=-180.0000; dip_dir=0;
slab2_cont='alu_slab2_dep_02.23.18_contours.in'; out11='Alu';
./extract_dip_east-west_2Contours $lon_min $lon_max $dip_depth $dip_depth_shall $dip_dir $smooth $box_width  $dip_int $slab2_cont $out11

lon_min=146.5; lon_max=155; dip_depth=-300.0000; dip_depth_shall=-260.0000; dip_dir=0;
slab2_cont='sol_slab2_dep_02.23.18_contours.in'; out12='Sol';
./extract_dip_east-west_2Contours $lon_min $lon_max $dip_depth $dip_depth_shall $dip_dir $smooth $box_width  $dip_int $slab2_cont $out12

text_dir='dip_catalogues/Slab2_two-depths'

if [ "$smooth" -eq 0 ]; then
	name="$text_dir/AllDips_2Contours.txt"
	name_full="$text_dir/AllDips_2Contours.full.txt"
else
	name="$text_dir/AllDips_2Contours.smoothed.txt"
	name_full="$text_dir/AllDips_2Contours.full.smoothed.txt"
fi

rm $text_dir/AllDips.txt 2> /dev/null
cat $text_dir/$out1.txt $text_dir/$out2.txt $text_dir/$out3.txt $text_dir/$out4.txt $text_dir/$out5.txt $text_dir/$out6.txt $text_dir/$out7.txt $text_dir/$out8.txt \
	$text_dir/$out9.txt $text_dir/$out10.txt $text_dir/$out11.txt $text_dir/$out12.txt > $name

rm $text_dir/AllDips.full.txt 2> /dev/null
cat $text_dir/$out1.full.txt $text_dir/$out2.full.txt $text_dir/$out3.full.txt $text_dir/$out4.full.txt $text_dir/$out5.full.txt $text_dir/$out6.full.txt $text_dir/$out7.full.txt \
	$text_dir/$out8.full.txt $text_dir/$out9.full.txt $text_dir/$out10.full.txt $text_dir/$out11.full.txt $text_dir/$out12.full.txt > $name_full

if [ "$filter" -eq 1 ]; then
	python remove_anomalous_regions.py $name $name_full
fi

rm $text_dir/$out1.full.txt $text_dir/$out2.full.txt $text_dir/$out3.full.txt $text_dir/$out4.full.txt $text_dir/$out5.full.txt $text_dir/$out6.full.txt $text_dir/$out7.full.txt \
	$text_dir/$out8.full.txt $text_dir/$out9.full.txt $text_dir/$out10.full.txt $text_dir/$out11.full.txt $text_dir/$out12.full.txt 2> /dev/null

rm $text_dir/$out1.txt $text_dir/$out2.txt $text_dir/$out3.txt $text_dir/$out4.txt $text_dir/$out5.txt $text_dir/$out6.txt $text_dir/$out7.txt $text_dir/$out8.txt \
	$text_dir/$out9.txt $text_dir/$out10.txt $text_dir/$out11.txt $text_dir/$out12.txt 2> /dev/null