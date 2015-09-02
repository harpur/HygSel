

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0006935","chemotaxis", 0.336,-0.116,-6.999, 5.224,-2.5342,0.830,0.000),
c("GO:0040011","locomotion", 1.437,-5.937,-0.652, 5.854,-2.5342,0.952,0.000),
c("GO:0048812","neuron projection morphogenesis", 0.032, 4.698,-0.850, 4.205,-3.5843,0.312,0.000),
c("GO:0051168","nuclear export", 0.078,-2.849, 2.955, 4.587,-1.7513,0.917,0.000),
c("GO:0071840","cellular component organization or biogenesis", 5.425,-4.716,-4.082, 6.431,-1.9496,0.954,0.000),
c("GO:0008202","steroid metabolic process", 0.040, 0.312, 6.985, 4.304,-1.7513,0.877,0.036),
c("GO:0010604","positive regulation of macromolecule metabolic process", 0.306,-5.495, 3.974, 5.182,-1.4540,0.619,0.037),
c("GO:0007033","vacuole organization", 0.009, 6.784, 0.152, 3.639,-1.7513,0.699,0.347),
c("GO:0006694","steroid biosynthetic process", 0.030, 1.196, 7.088, 4.178,-1.7513,0.859,0.423),
c("GO:0030030","cell projection organization", 0.378, 5.739, 1.058, 5.274,-2.3369,0.626,0.443),
c("GO:0009605","response to external stimulus", 1.384, 0.921,-7.190, 5.838,-1.9047,0.881,0.444),
c("GO:0007040","lysosome organization", 0.003, 6.579, 1.156, 3.214,-1.7513,0.711,0.471),
c("GO:0022412","cellular process involved in reproduction in multicellular organism", 0.024, 0.931,-2.516, 4.076,-1.7513,0.603,0.586),
c("GO:0016043","cellular component organization", 4.287, 6.147,-0.236, 6.329,-1.9955,0.627,0.632));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);



# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("C:/Users/path_to_your_file/revigo-plot.pdf");
