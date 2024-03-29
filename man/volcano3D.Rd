% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/volcano3D.R
\name{volcano3D}
\alias{volcano3D}
\title{Three-Dimensional Volcano Plot}
\usage{
volcano3D(
  polar,
  type = 1,
  label_rows = c(),
  label_size = 14,
  arrow_length = 100,
  colour_code_labels = FALSE,
  label_colour = "black",
  grid_colour = "grey80",
  grid_width = 2,
  grid_options = NULL,
  axis_colour = "black",
  axis_width = 2,
  marker_size = 3,
  marker_outline_width = 0,
  marker_outline_colour = "white",
  z_axis_title_offset = 1.2,
  z_axis_title_size = 12,
  z_axis_angle = 0.5,
  radial_axis_title_size = 14,
  radial_axis_title_offset = 1.2,
  xy_aspectratio = 1,
  z_aspectratio = 0.8,
  camera_eye = list(x = 0.9, y = 0.9, z = 0.9),
  ...
)
}
\arguments{
\item{polar}{Object of S4 class 'volc3d' following call to either
\code{\link{polar_coords}}, \code{\link{deseq_polar}} or
\code{\link{voom_polar}}}

\item{type}{Either \code{1} or \code{2} specifying type of polar coordinates: \code{1} =
Z-scaled, \code{2} = unscaled (equivalent to log2 fold change for gene
expression).}

\item{label_rows}{A vector of row names or numbers to label}

\item{label_size}{font size for labels (default 14).}

\item{arrow_length}{The length of label arrows (default 100)}

\item{colour_code_labels}{Logical whether label annotations should be colour
coded. If \code{FALSE} \code{label_colour} is used.}

\item{label_colour}{HTML colour of annotation labels if not colour coded.}

\item{grid_colour}{The colour of the cylindrical grid (default "grey80")}

\item{grid_width}{The width of the grid lines (default 2)}

\item{grid_options}{Optional list of additional arguments to pass to
\code{\link{polar_grid}}, eg. \code{z_axis_ticks} and \code{r_axis_ticks}}

\item{axis_colour}{The colour of the grid axes and labels (default "black")}

\item{axis_width}{The width of axis lines (default 2)}

\item{marker_size}{Size of the markers (default 3)}

\item{marker_outline_width}{Width for marker outline (default 0 means no
outline)}

\item{marker_outline_colour}{Colour for marker outline (default white)}

\item{z_axis_title_offset}{The position scaling between grid and z axis title
(default=1.2)}

\item{z_axis_title_size}{The font size for the z axis title (default=12)}

\item{z_axis_angle}{Angle in radians for the position of z axis (default
0.5)}

\item{radial_axis_title_size}{The font size for the radial (default=15)}

\item{radial_axis_title_offset}{The position scaling between grid and radial
axis title (default=1.2)}

\item{xy_aspectratio}{The aspect ratio for the xy axis compared to z (default
1). Increasing this makes the grid wider in the plot window.}

\item{z_aspectratio}{The aspect ratio for the z axis compared to x and y
(default 0.8). Decreasing this makes the plot appear more squat.}

\item{camera_eye}{The (x,y,z) components of the start 'eye' camera vector.
This vector determines the view point about the origin of this scene.}

\item{...}{Optional arguments passed to \code{\link[plotly:plot_ly]{plot_ly}}}
}
\value{
Returns a cylindrical 3D plotly plot featuring variables on a
tri-axis radial graph with the -log10(multi-group test p-value) on the
z-axis
}
\description{
Plots the three-way comparisons of variables such as gene expression data in
3D space using plotly. x, y position represents polar position on 3 axes
representing the amount each variable or gene tends to each of the 3
categories. The z axis represents -log10 P value for the one-way test
comparing each variable across the 3 groups.
}
\examples{
data(example_data)
syn_polar <- polar_coords(outcome = syn_example_meta$Pathotype,
                          data = t(syn_example_rld))
volcano3D(syn_polar)

}
\references{
Lewis, Myles J., et al. (2019).
\href{https://pubmed.ncbi.nlm.nih.gov/31461658/}{
Molecular portraits of early rheumatoid arthritis identify clinical and
treatment response phenotypes.}
\emph{Cell reports}, \strong{28}:9
}
\seealso{
\code{\link{polar_coords}}
}
