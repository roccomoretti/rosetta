# -*- tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
# vi: set ts=2 noet:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

check_setup()
feature_analyses <- c(feature_analyses, new("FeaturesAnalysis",
id = "cation_pi_eqpoldens_scOrbCpi_Hpol_chiBDHO_cosDHO",
author = "Matthew O'Meara, Steven Combs",
brief_description = "",
feature_reporter_dependencies = c("OrbitalFeatures"),
run=function(self, sample_sources, output_dir, output_formats){


						
#################################################################################
###########################All cation-pi at 2.0A#################################
#################################################################################
sele <- "
SELECT
  resNum1,
  resName1,
  resNum2,
  resName2,
  OrbName1,
  cosDHO,
  OrbHdist,
  chiBDHO,
  htype2
FROM
  Hpol_orbital
WHERE
  (((resName2 = 'LYS' OR resName2 = 'ARG' ) AND (resName1 = 'PHE' OR resName1 = 'TRP' OR resName1 = 'TYR')))  AND 
  OrbHdist < 2.0 AND 
  OrbName1 = 'C.pi.sp2' AND
  htype2 = 'Hpol' AND 
  ABS(resNum1 - resNum2) > 5;"
f <- query_sample_sources(sample_sources, sele)
f <- transform(f,
  capx = 2*sin(acos(cosDHO)/2)*cos(chiBDHO),
  capy = 2*sin(acos(cosDHO)/2)*sin(chiBDHO))

capx_limits <- c(-1.5,1.5)
capy_limits <- capx_limits
plot_id = "cation_pi_eqpoldens_scOrbCpi_Hpol_chiBDHO_cosDHO_all_2A"
f_first <- f[ f$sample_source == levels(sample_sources$sample_source), ]

#ggplot(data=f_first) + theme_bw() +
#  theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
#  stat_density2d(
#    aes(x=capx, y=capy, fill=..density..), geom="tile", contour=FALSE ) +
#  polar_equal_area_grids_bw() +
#  ggtitle(
#    paste("Cation pi chiBDHO vs DHO Angles with Sequence Separation > 5\n",
#    "Sidechain Donors to Sidechain sp2 Acceptors, Equal Coordinate Projection\n",
#    "LYS+ARG to PHE+TYR+TRP at 2.0", sep="")) +
#  scale_x_continuous(
#    '2*sin(DHO/2) * cos(chiBDHO)', limits=capx_limits, breaks=c(-1, 0, 1)) +
#  scale_y_continuous(
#    '2*sin(DHO/2) * sin(chiBDHO)', limits=capy_limits, breaks=c(-1, 0, 1)) +
#  coord_fixed(ratio = 1) +
#  scale_fill_gradientn('Density', colour=jet.colors(10))
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)


ggplot(data=f_first) + theme_bw() +
		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
		stat_density2d(
				aes(x=capx, y=capy, fill=..density..), geom="tile", contour=FALSE ) +
		polar_equal_area_grids_bw() +
		ggtitle(
						paste("Cation pi chiBDHO vs DHO Angles with Sequence Separation > 5\n",
    "Sidechain Donors to Sidechain sp2 Acceptors, Equal Coordinate Projection\n",
    "LYS+ARG to PHE+TYR+TRP at 2.0)", sep="")) +
		scale_x_continuous(
				'2*sin(DHO/2) * cos(chiBDHO)', limits=capx_limits, breaks=c(-1, 0, 1)) +
		scale_y_continuous(
				'2*sin(DHO/2) * sin(chiBDHO)', limits=capy_limits, breaks=c(-1, 0, 1)) +
		coord_fixed(ratio = 1) +
		scale_fill_gradientn('Density', colours=jet.colors(10))
save_plots(self, plot_id, sample_sources, output_dir, output_formats)



})) # end FeaturesAnalysis

#################################################################################
###########################Phe-lys cation-pi at 2.0A#############################
#################################################################################
#sele <- "
#		SELECT
#		resNum1,
#		resName1,
#		resNum2,
#		resName2,
#		OrbName1,
#		cosDHO,
#		OrbHdist,
#		chiBDHO,
#		htype2
#		FROM
#		Hpol_orbital
#		WHERE
#		(((resName2 = 'LYS' ) AND (resName1 = 'PHE' )))  AND
#		OrbHdist < 2 AND
#		OrbName1 = 'C.pi.sp2' AND
#		htype2 = 'Hpol' AND
#		ABS(resNum1 - resNum2) > 5;"
#
#f <- query_sample_sources(sample_sources, sele)
#
#f <- transform(f,
#		capx = 2*sin(acos(cosDHO)/2)*cos(chiBDHO),
#		capy = 2*sin(acos(cosDHO)/2)*sin(chiBDHO))
#
#capx_limits <- c(-1.5,1.5)
#capy_limits <- capx_limits
#
#
#plot_id = "chiBDHO_cosDHO_eqpoldens_cation_pi_cpisp2_lys_phe_2A"
#
#f_first <- f[ f$sample_source == levels(sample_sources$sample_source), ]
#
#ggplot(data=f_first) + theme_bw() +
#		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
#		stat_density2d(
#				aes(x=capx, y=capy, fill=..density..), geom="tile", contour=FALSE ) +
#		polar_equal_area_grids_bw() +
#		ggtitle(
#						paste("Cation pi chiBDHO vs DHO Angles with Sequence Separation > 5\n",
#								"Sidechain Donors to Sidechain sp2 Acceptors, Equal Coordinate Projection\n",
#								"LYS to PHE at 2.0", sep="")) +
#		scale_x_continuous(
#				'2*sin(DHO/2) * cos(chiBDHO)', limits=capx_limits, breaks=c(-1, 0, 1)) +
#		scale_y_continuous(
#				'2*sin(DHO/2) * sin(chiBDHO)', limits=capy_limits, breaks=c(-1, 0, 1)) +
#		coord_fixed(ratio = 1) +
#		scale_fill_gradientn('Density', colour=jet.colors(10))
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
#
##################################################################################
############################Phe-arg cation-pi at 2.0A#############################
##################################################################################
#sele <- "
#		SELECT
#		resNum1,
#		resName1,
#		resNum2,
#		resName2,
#		OrbName1,
#		cosDHO,
#		OrbHdist,
#		chiBDHO,
#		htype2
#		FROM
#		Hpol_orbital
#		WHERE
#		(((resName2 = 'ARG' ) AND (resName1 = 'PHE' )))  AND
#		OrbHdist < 2 AND
#		OrbName1 = 'C.pi.sp2' AND
#		htype2 = 'Hpol' AND
#		ABS(resNum1 - resNum2) > 5;"
#
#f <- query_sample_sources(sample_sources, sele)
#
#f <- transform(f,
#		capx = 2*sin(acos(cosDHO)/2)*cos(chiBDHO),
#		capy = 2*sin(acos(cosDHO)/2)*sin(chiBDHO))
#
#capx_limits <- c(-1.5,1.5)
#capy_limits <- capx_limits
#
#
#plot_id = "chiBDHO_cosDHO_eqpoldens_cation_pi_cpisp2_arg_phe_2A"
#
#f_first <- f[ f$sample_source == levels(sample_sources$sample_source), ]
#
#ggplot(data=f_first) + theme_bw() +
#		theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
#		stat_density2d(
#				aes(x=capx, y=capy, fill=..density..), geom="tile", contour=FALSE ) +
#		polar_equal_area_grids_bw() +
#		ggtitle(
#						paste("Cation pi chiBDHO vs DHO Angles with Sequence Separation > 5\n",
#								"Sidechain Donors to Sidechain sp2 Acceptors, Equal Coordinate Projection\n",
#								"ARG to PHE at 2.0", sep="")) +
#		scale_x_continuous(
#				'2*sin(DHO/2) * cos(chiBDHO)', limits=capx_limits, breaks=c(-1, 0, 1)) +
#		scale_y_continuous(
#				'2*sin(DHO/2) * sin(chiBDHO)', limits=capy_limits, breaks=c(-1, 0, 1)) +
#		coord_fixed(ratio = 1) +
#		scale_fill_gradientn('Density', colour=jet.colors(10))
#save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
#
#
####Phe-Lys 3.5
##sele <- "
##SELECT
##  resNum1,
##  resName1,
##  resNum2,
##  resName2,
##  OrbName1,
##  cosDHO,
##  OrbHdist,
##  chiBDHO,
##  htype2
##FROM
##  Hpol_orbital
##WHERE
##  (((resName2 = 'LYS' ) AND (resName1 = 'PHE' )))  AND
##  OrbHdist < 3.5 AND
##  OrbName1 = 'C.pi.sp2' AND
##  htype2 = 'Hpol' AND
##  ABS(resNum1 - resNum2) > 5;"
##
##f <- query_sample_sources(sample_sources, sele)
##
##f <- transform(f,
##  capx = 2*sin(acos(cosDHO)/2)*cos(chiBDHO),
##  capy = 2*sin(acos(cosDHO)/2)*sin(chiBDHO))
##
##capx_limits <- c(-1.5,1.5)
##capy_limits <- capx_limits
##
##
##plot_id = "chiBDHO_cosDHO_eqpoldens_cation_pi_cpisp2_lys_phe_3.5A"
##
##f_first <- f[ f$sample_source == levels(sample_sources$sample_source), ]
##
##ggplot(data=f_first) + theme_bw() +
##  theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
##  stat_density2d(
##    aes(x=capx, y=capy, fill=..density..), geom="tile", contour=FALSE ) +
##  polar_equal_area_grids_bw() +
##  ggtitle(
##    paste("Cation pi chiBDHO vs DHO Angles with Sequence Separation > 5\n",
##    "Sidechain Donors to Sidechain sp2 Acceptors, Equal Coordinate Projection\n",
##    "LYS to PHE at 3.5", sep="")) +
##  scale_x_continuous(
##    '2*sin(DHO/2) * cos(chiBDHO)', limits=capx_limits, breaks=c(-1, 0, 1)) +
##  scale_y_continuous(
##    '2*sin(DHO/2) * sin(chiBDHO)', limits=capy_limits, breaks=c(-1, 0, 1)) +
##  coord_fixed(ratio = 1) +
##  scale_fill_gradientn('Density', colour=jet.colors(10))
##save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
####phe - arg 3.5
##sele <- "
##SELECT
##  resNum1,
##  resName1,
##  resNum2,
##  resName2,
##  OrbName1,
##  cosDHO,
##  OrbHdist,
##  chiBDHO,
##  htype2
##FROM
##  Hpol_orbital
##WHERE
##  (((resName2 = 'ARG' ) AND (resName1 = 'PHE' )))  AND
##  OrbHdist < 3.5 AND
##  OrbName1 = 'C.pi.sp2' AND
##  htype2 = 'Hpol' AND
##  ABS(resNum1 - resNum2) > 5;"
##
##f <- query_sample_sources(sample_sources, sele)
##
##f <- transform(f,
##  capx = 2*sin(acos(cosDHO)/2)*cos(chiBDHO),
##  capy = 2*sin(acos(cosDHO)/2)*sin(chiBDHO))
##
##capx_limits <- c(-1.5,1.5)
##capy_limits <- capx_limits
##
##
##plot_id = "chiBDHO_cosDHO_eqpoldens_cation_pi_cpisp2_arg_phe_3.5A"
##
##f_first <- f[ f$sample_source == levels(sample_sources$sample_source), ]
##
##ggplot(data=f_first) + theme_bw() +
##  theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
##  stat_density2d(
##    aes(x=capx, y=capy, fill=..density..), geom="tile", contour=FALSE ) +
##  polar_equal_area_grids_bw() +
##  ggtitle(
##    paste("Cation pi chiBDHO vs DHO Angles with Sequence Separation > 5\n",
##    "Sidechain Donors to Sidechain sp2 Acceptors, Equal Coordinate Projection\n",
##    "ARG to PHE at 3.5", sep="")) +
##  scale_x_continuous(
##    '2*sin(DHO/2) * cos(chiBDHO)', limits=capx_limits, breaks=c(-1, 0, 1)) +
##  scale_y_continuous(
##    '2*sin(DHO/2) * sin(chiBDHO)', limits=capy_limits, breaks=c(-1, 0, 1)) +
##  coord_fixed(ratio = 1) +
##  scale_fill_gradientn('Density', colour=jet.colors(10))
##save_plots(self, plot_id, sample_sources, output_dir, output_formats)
#
######cation pi interactions all at 3.5 A
##sele <- "
##SELECT
##  resNum1,
##  resName1,
##  resNum2,
##  resName2,
##	OrbName1,
##	cosDHO,
##  OrbHdist,
##	chiBDHO,
##	htype2
##FROM
##  Hpol_orbital
##WHERE
##  (((resName2 = 'LYS' OR resName2 = 'ARG' ) AND (resName1 = 'PHE' OR resName1 = 'TRP' OR resName1 = 'TYR')))  AND 
##	OrbHdist < 3.5 AND 
##	OrbName1 = 'C.pi.sp2' AND
##	htype2 = 'Hpol' AND 
##	ABS(resNum1 - resNum2) > 5;"
##f <- query_sample_sources(sample_sources, sele)
##f <- transform(f,
##  capx = 2*sin(acos(cosDHO)/2)*cos(chiBDHO),
##  capy = 2*sin(acos(cosDHO)/2)*sin(chiBDHO))
##
##capx_limits <- c(-1.5,1.5)
##capy_limits <- capx_limits
##plot_id = "chiBDHO_cosDHO_eqpoldens_cation_pi_cpisp2_all_interactions_3.5A"
##f_first <- f[ f$sample_source == levels(sample_sources$sample_source), ]
##
##ggplot(data=f_first) + theme_bw() +
##  theme(panel.background=element_rect(fill="#00007F", colour="#00007F")) +
##  stat_density2d(
##    aes(x=capx, y=capy, fill=..density..), geom="tile", contour=FALSE ) +
##  polar_equal_area_grids_bw() +
##  ggtitle(
##    paste("Cation pi chiBDHO vs DHO Angles with Sequence Separation > 5 && <3.5A\n",
##    "Sidechain Donors to Sidechain sp2 Acceptors, Equal Coordinate Projection\n",
##		"LYS+ARG to PHE+TYR+TRP", sep="")) +
##  scale_x_continuous(
##    '2*sin(DHO/2) * cos(chiBDHO)', limits=capx_limits, breaks=c(-1, 0, 1)) +
##  scale_y_continuous(
##    '2*sin(DHO/2) * sin(chiBDHO)', limits=capy_limits, breaks=c(-1, 0, 1)) +
##  coord_fixed(ratio = 1) +
##  scale_fill_gradientn('Density', colour=jet.colors(10))
##save_plots(self, plot_id, sample_sources, output_dir, output_formats)

