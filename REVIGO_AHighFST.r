

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
revigo.data <- rbind(c("GO:0000003","reproduction", 0.214,-4.872, 4.273, 5.028,-1.9430,1.000,0.000),
c("GO:0007610","behavior", 0.052,-1.630, 7.317, 4.416,-2.6278,0.887,0.000),
c("GO:0022414","reproductive process", 0.121, 2.370,-5.480, 4.781,-1.9430,0.955,0.000),
c("GO:0023052","signaling", 3.838,-2.180,-4.089, 6.281,-1.8533,0.957,0.000),
c("GO:0032501","multicellular organismal process", 0.790,-5.030, 1.781, 5.594,-2.4086,0.955,0.000),
c("GO:0032502","developmental process", 1.387,-4.634,-1.193, 5.839,-1.4920,0.956,0.000),
c("GO:0048468","cell development", 0.132, 7.416, 1.235, 4.818,-3.2040,0.389,0.000),
c("GO:0048610","cellular process involved in reproduction", 0.159,-2.327,-0.885, 4.849,-2.3523,0.955,0.000),
c("GO:0050896","response to stimulus", 8.818, 0.248,-4.213, 6.642,-1.8008,0.959,0.000),
c("GO:0007154","cell communication", 4.358, 5.276,-5.124, 6.336,-1.8533,0.891,0.103),
c("GO:0007626","locomotory behavior", 0.015, 0.073, 8.049, 3.865,-2.1044,0.825,0.270),
c("GO:0009605","response to external stimulus", 1.384,-0.987, 7.830, 5.838,-1.3718,0.876,0.374),
c("GO:0007476","imaginal disc-derived wing morphogenesis", 0.002, 5.538, 3.098, 2.960,-2.0892,0.281,0.495),
c("GO:0007552","metamorphosis", 0.003, 6.976, 0.497, 3.134,-1.5253,0.499,0.507),
c("GO:0007417","central nervous system development", 0.052, 5.921, 3.368, 4.417,-2.5690,0.277,0.612),
c("GO:0003006","developmental process involved in reproduction", 0.081, 7.124, 1.006, 4.603,-2.3523,0.349,0.614),
c("GO:0032504","multicellular organism reproduction", 0.052, 5.272, 5.655, 4.416,-1.9430,0.357,0.636),
c("GO:0048736","appendage development", 0.015, 5.918, 2.833, 3.870,-2.0892,0.335,0.655),
c("GO:0016319","mushroom body development", 0.000, 6.157, 3.740, 2.391,-2.4019,0.378,0.679),
c("GO:0008345","larval locomotory behavior", 0.000, 2.914, 5.950, 1.924,-1.9975,0.446,0.682),
c("GO:0048569","post-embryonic organ development", 0.015, 6.263, 3.216, 3.880,-1.8389,0.317,0.688),
c("GO:0030537","larval behavior", 0.000, 3.135, 6.154, 2.061,-1.9975,0.440,0.692));

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
