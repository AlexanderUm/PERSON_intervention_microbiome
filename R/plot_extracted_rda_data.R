#-------------------------------------------------------------------------------
# Plot dbRDA scatterpots 
#-------------------------------------------------------------------------------

plot_extracted_rda_data <- function(extracted_data, 
                                    connect_observations = TRUE, 
                                    add_elepses = TRUE, 
                                    group_col,
                                    observ_time_group_col, 
                                    time_col, 
                                    color_vec = NULL, 
                                    shape_vec = NULL, 
                                    stat_text = NULL,
                                    dists_named_vec = NULL) {
  
  # Plot individual plots 
  ind.plots <- list()
  
  for(i in names(extracted_data)) {
    
    # Plot title 
    if(is.null(dists_named_vec)) {
      
      p.title <- paste0("Distance: ", i)
      
    } else {
      
      p.title <- paste0("Distance: ", 
                        names(dists_named_vec)[dists_named_vec == i])
      
    }
      
    # extract data 
    p.df.inst <- extracted_data[[i]]$main %>% 
                    arrange(.data[[observ_time_group_col]], 
                            .data[[time_col]])
    
    # Base plot
    p.out <- ggplot(p.df.inst,
                     aes(x = dbRDA1, 
                         y = dbRDA2)) +
              ggtitle(p.title) +
              coord_fixed() +
              theme_bw() + 
              theme(axis.line = element_line(color='black'),
                    plot.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) + 
              xlab(extracted_data[[i]]$var_expl[1]) +
              ylab(extracted_data[[i]]$var_expl[2]) +
              theme(plot.title = element_text(size=11, 
                                              face="italic"))
    
      # Add ellipse layer
      if(add_elepses) {
      
        p.out  <- p.out +
                    stat_ellipse(aes(fill = .data[[group_col]]), 
                                 geom = "polygon", 
                                 alpha = 0.15, 
                                 level = 0.90, 
                                 size =2) 
      }
    
      # Add path layer
      if(connect_observations) {
      
        p.out <- p.out + 
                geom_path(aes(group = .data[[observ_time_group_col]]),
                          alpha = 1, 
                          color = "gray",
                          linewidth = 0.5) 
      }
    
      
      # Add custom color
      if(!is.null(color_vec)) {
        
        p.out <- p.out + 
                  scale_color_manual(values = color_vec) +
                  scale_fill_manual(values = color_vec)
      }
      
      # Add custom shape
      if(!is.null(shape_vec)) {
        p.out <- p.out + 
                  scale_shape_manual(values = shape_vec)
      }
    
      # Add statistics text
      if(!is.null(stat_text)) {
        
        p.text.inst <- stat_text %>% 
                        filter(Distance == i) %>% 
                        mutate(x = min(p.df.inst$dbRDA1), 
                               y = max(p.df.inst$dbRDA2) + 
                                   (max(p.df.inst$dbRDA2) - 
                                      min(p.df.inst$dbRDA2))*0.1)
        
        
        p.out <- p.out + 
                  geom_text(data = p.text.inst, 
                            aes(x = x, 
                                y = y, 
                                label = Text), 
                            hjust = 0, parse = T, 
                            size = 3.5)
        
      }
      
      # Add points layer 
      p.out <- p.out + 
                geom_point(aes(color = .data[[group_col]],
                               shape = .data[[time_col]]),
                           alpha = 1, 
                           size = 2)
      
      ind.plots[[i]] <- p.out
      
  }
  
  #-------------------------------------------------------------------------------
  # Combine plots 
  #-------------------------------------------------------------------------------
  ind.plots.f <- lapply(ind.plots,
                        function(x){x+theme(legend.position="none")})
  
  legd.inst <- get_legend(ind.plots[[1]])
  
  p.grid.1 <- plot_grid(plotlist = ind.plots.f, 
                        ncol = 2)
  
  p.grid.2 <- plot_grid(p.grid.1, 
                        legd.inst, 
                        ncol = 2, 
                        rel_widths = c(0.85, .15))
  
  res.out <- list("Comb" = p.grid.2, 
                  "Ind" = ind.plots, 
                  "Legend" = legd.inst)
  
  return(res.out)
  
}