library(GGally)


my_rosenbrock<-function(x,a=0,b=100,c=0.1,log=FALSE){
  
  d<-length(x)
  y<-x^2
  
  term1<-b*sum((x[2:d]-y[1:(d-1)])^2)
  term2<-sum((a-x[1:(d-1)])^2)
  
  res<- -0.5/c*(term1+term2)
  
  if(log==TRUE)
    return(res)
  else
    return(exp(res))
}




a <-0
b <-100
c <- 0.1

n<-1e+6
d<-3

x<-matrix(rnorm(n*d),ncol=d)

x[,1]<-x[,1]*sqrt(c)+a
for(i in 2:max(d-1,2)){
  x[,i]<-x[,i]*sqrt(c/(b+1))+(b*x[,i-1]^2+a)/(b+1)
}
if(d>2)
  x[,d]<-x[,d]*sqrt(c/b)+x[,d-1]^2


# OLD PLOT
# # Modify the column names to use expression with subscripts
# colnames(x) <- c(expression(x[1]), expression(x[2]), expression(x[3]))
# 
# # Define upper panel with black contour lines
# upperfun <- function(data, mapping) {
#   ggplot(data = data, mapping = mapping) +
#     geom_density2d(color = "black") +
#     theme_bw() +
#     theme(
#       axis.text = element_text(size = 8), 
#       axis.ticks = element_line()
#     )
# }
# 
# # Define lower panel with black points
# lowerfun <- function(data, mapping) {
#   ggplot(data = data, mapping = mapping) +
#     geom_point(color = "black", alpha = 0.5) +
#     theme_bw() +
#     theme(
#       axis.text = element_text(size = 8), 
#       axis.ticks = element_line()
#     )
# }
# 
# # Diagonal panel: Black density curves with white fill
# diagfun <- function(data, mapping) {
#   ggplot(data = data, mapping = mapping) +
#     geom_density(fill = "white", color = "black") +
#     theme_bw() +
#     theme(
#       axis.text = element_text(size = 8), 
#       axis.ticks = element_line()
#     )
# }
# 
# # Generate the ggpairs plot with BW styling and LaTeX-style labels
# ggpairs(as.data.frame(x),
#         upper = list(continuous = wrap(upperfun)),
#         lower = list(continuous = wrap(lowerfun)),
#         diag = list(continuous = wrap(diagfun)),
#         labeller = label_parsed)
# 
# ggsave("density_ex2_bw.png", width = 10, height = 8, dpi = 300, units = "in")
