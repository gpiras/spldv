\name{NEWS}
\title{spldv News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 0.2.2}{
  \subsection{BUG FIXES AND IMPROVEMENTS}{
    \itemize{
      \item Fixed a bug in the computation of the spatial multiplier matrix used by marginal effects.
      \item For continuous variables, \code{impacts(..., result = "from.region", dydx = "exact")} and \code{impacts(..., result = "cumulative", dydx = "exact")} now use analytical derivatives rather than finite differences. This improves numerical stability and makes the \code{"exact"} option consistent across average, from-region, and cumulative impacts.
      \item \code{impacts()} now stops with an informative error when the fitted formula contains transformations that are not currently supported for spatial marginal effects, including \code{poly()}, \code{scale()}, \code{bs()}, and \code{ns()}. For quadratic effects, users should write explicit terms such as \code{x + I(x^2)}. For scaled variables, users should create the transformed variable in the data before fitting the model.
      \item Discrete variables, including factors, ordered factors, and logical variables, continue to be evaluated as predictive contrasts rather than continuous derivatives. This includes formulas that use inline factors such as \code{factor(x)} and interactions involving such factors.
    }
  }
  \subsection{DOCUMENTATION}{
    \itemize{
      \item Expanded the help page for \code{impacts()} and its \code{summary()} method. The documentation now explains average direct, indirect, and total impacts; from-region impacts; cumulative impacts; the role of \code{from.unit}; the interpretation of \code{Q}; the difference between \code{dydx = "exact"} and \code{dydx = "numeric"}; and the distinction between continuous derivatives and discrete contrasts.
      \item Added examples illustrating average impacts, from-region impacts, cumulative impacts, Monte Carlo and Delta-method standard errors, finite changes through \code{change = c(lower, upper)}, and impacts for models containing inline factor terms.
    }
  }
  \subsection{TESTS AND MAINTENANCE}{
    \itemize{
      \item Added regression tests for inline factor terms in average, from-region, and cumulative impacts.
      \item Ignored local diagnostic and check artifacts so they are not included in package builds.
    }
  }
}

\section{CHANGES IN VERSION 0.2.0}{
\subsection{NEW FEATURES}{
    \itemize{
      \item Includes prediction function for GMM, LGMM, and RIS estimators.
      \item Includes the argument fastmom argument in sbinaryGMM function which speed up the computation time.
      \item New functionalities for computing marginal effects. 
      \item Fixed a bug in the gradient of GMM models. 
      \item atmeans argument for marginal effects is deprecated.
      \item Fixed overflow in the S matrix.
      \item Improved documentation.
    }
  }  
}

\section{CHANGES IN VERSION 0.1.4}{
\subsection{NEW FEATURES}{
    \itemize{
      \item Includes RIS estimator with numerical gradient and marginal effects.
    }
  }  
}

\section{CHANGES IN VERSION 0.1.3}{
\subsection{NEW FEATURES}{
    \itemize{
      \item Includes citation for JSS paper.
    }
  }  
}

\section{CHANGES IN VERSION 0.1.2}{
\subsection{NEW FEATURES}{
    \itemize{
      \item Improved documentation.
      \item Fixed bug with impacts.
    }
  }  
}

\section{CHANGES IN VERSION 0.1.1}{
\subsection{NEW FEATURES}{
    \itemize{
      \item Fixed bug with constrained optimization (JSS). 
      \item Fixed minor bugs with controls for maxLik.
      \item Improved documentation. 
      \item Added impacts() function for marginal effects. 
      \item Added tol.solve for solve function. 
    }
  }  
}
