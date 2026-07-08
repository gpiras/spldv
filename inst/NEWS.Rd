\name{NEWS}
\title{spldv News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 0.2.2}{
\subsection{NEW FEATURES}{
    \itemize{
      \item Fixes bug in S matrix.
      \item \code{impacts()} now stops with an informative error when the model formula contains transformed terms currently unsupported for spatial marginal effects, such as \code{poly()}, \code{scale()}, \code{bs()}, or \code{ns()}. For quadratic effects, users should write explicit terms such as \code{x + I(x^2)}. For scaled variables, users should create the transformed variable in the data before fitting the model.
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

