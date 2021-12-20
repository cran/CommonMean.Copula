#' Estimate bivariate common mean vector under copula models
#'
#' @param Y1 Outcome 1
#' @param Y2 Outcome 2
#' @param Sigma1 Standard deviation of outcome 1.
#' @param Sigma2 Standard deviation of outcome 2.
#' @param rho Correlation coefficient between outcomes.
#' @param copula The copula to be used with possible options \code{"Clayton"}, \code{"Gumbel"}, \code{"Frank"}, \code{"FGM"}, and \code{"normal"}.
#' @description Estimate the common mean vector under copula models with known correlation. A maximum likelihood estimation procedure is employed. See Shih et al. (2019) and Shih et al. (2021) for details under the Farlie-Gumbel-Morgenstern (FGM) and general copulas, respectively.
#' @details We apply \code{"optim"} routine to maximize the log-likelihood function. In addition, boundary corrected correlations will be used (Shih et al., 2019).
#' @return \item{Outcome 1}{Outcome 1.}
#' \item{Outcome 2}{Outcome 2.}
#' \item{Correlation}{Correlation coefficient between outcomes.}
#' \item{Sample size}{Sample size.}
#' \item{Copula}{Selected copula.}
#' \item{Copula parameter}{Copula parameter.}
#' \item{Corrected correlation}{Boundary corrected correlations.}
#' \item{CommonMean 1}{Estimation results of outcome 1.}
#' \item{CommonMean 2}{Estimation results of outcome 2.}
#' \item{V}{Covariance matrix of the common mean vector estimate.}
#' \item{Log-likelihood values}{Fitted log-likelihood values.}
#'
#' @note When \code{rho} is 1 or -1, there are some computational issues since the copula parameter may correspond to infinite or negative infinite under some copulas. For the Clayton copula, if \code{rho} > 0.95, it will be approximated by 0.95. For the Frank copula, if \code{rho} > 0.95 or \code{rho} < -0.95, it will be approximated by 0.95 or -0.95, respectively.
#' @references Shih J-H, Konno Y, Chang Y-T, Emura T (2019) Estimation of a common mean vector in bivariate meta-analysis under the FGM copula, Statistics 53(3): 673-95.
#' @references Shih J-H, Konno Y, Emura T (2021-) Copula-based estimation methods for a common mean vector for bivariate meta-analyses, under review.
#' @importFrom stats optim pnorm qnorm uniroot
#' @importFrom pracma integral2
#' @importFrom mvtnorm dmvnorm
#' @export
#'
#' @examples
#' library(CommonMean.Copula)
#' Y1 = c(35,25,30,50,60) # outcome 1
#' Y2 = c(30,30,50,65,40) # outcome 2
#' Sigma1 = c(1.3,1.4,1.5,2.0,1.8) # SE of outcome 1
#' Sigma2 = c(1.7,1.9,2.5,2.2,1.8) # SE of outcome 2
#' rho = c(0.4,0.7,0.6,0.7,0.6) # correlation between two outcomes
#' CommonMean.Copula(Y1,Y2,Sigma1,Sigma2,rho) # input

CommonMean.Copula = function(Y1,Y2,Sigma1,Sigma2,rho,copula = "Clayton") {

  n = length(Y1)

  if (copula == "Clayton") {

    c_density = function(u,v,theta) {

      (theta+1)*(u*v)^(-theta-1)*(u^-theta+v^-theta-1)^(-1/theta-2)

    }

    ### boundary correction
    rho.BC = pmax(0.001,rho)
    rho.BC.temp = pmin(0.95,rho.BC)

    uni_cor = function(rho) {

      inv_f = function(theta) {

        cor_f = function(u,v) {qnorm(u)*qnorm(v)*c_density(u,v,theta)}
        integral2(cor_f,0,1,0,1)$Q-rho

      }
      uniroot(inv_f,c(0.001,12),extendInt = "yes")$root

    }
    theta.vec = sapply(rho.BC.temp,uni_cor)

  }
  if (copula == "Gumbel") {

    c_density = function(u,v,theta) {

      a1 = 1/u*(-log(u))^(theta-1)
      a2 = 1/v*(-log(v))^(theta-1)
      a3 = exp(-((-log(u))^theta+(-log(v))^theta)^(1/theta))
      a4 = ((-log(u))^theta+(-log(v))^theta)^(1/theta-2)
      a5 = (theta-1)+((-log(u))^theta+(-log(v))^theta)^(1/theta)

      a1*a2*a3*a4*a5

    }

    ### boundary correction
    rho.BC = pmax(0.001,rho)

    uni_cor = function(rho) {

      inv_f = function(theta) {

        cor_f = function(u,v) {qnorm(u)*qnorm(v)*c_density(u,v,theta)}
        integral2(cor_f,0,1,0,1)$Q-rho

      }
      uniroot(inv_f,c(1.001,6),extendInt = "yes")$root

    }
    theta.vec = sapply(rho.BC,uni_cor)

  }
  if (copula == "Frank") {

    c_density = function(u,v,theta) {

      theta*(1-exp(-theta))*exp(-theta*(u+v))/((1-exp(-theta))-(1-exp(-theta*u))*(1-exp(-theta*v)))^2

    }

    ### boundary correction
    rho.BC = rho
    rho.BC.temp = pmin(0.95,pmax(-0.95,rho.BC))

    uni_cor = function(rho) {

      inv_f = function(theta) {

        cor_f = function(u,v) {qnorm(u)*qnorm(v)*c_density(u,v,theta)}
        integral2(cor_f,0,1,0,1)$Q-rho

      }
      uniroot(inv_f,c(-30,30),extendInt = "yes")$root

    }
    theta.vec = sapply(rho.BC.temp,uni_cor)

  }
  if (copula == "FGM") {

    c_density = function(u,v,theta) {

      1+theta*(1-2*u)*(1-2*v)

    }
    theta.vec = pmin(1,pmax(rho*pi,-1))
    rho.BC = theta.vec/pi

  }

  general_LLH = function(para) {

    Mu1 = para[1]
    Mu2 = para[2]

    Phi1 = pnorm((Y1-Mu1)/Sigma1)
    Phi2 = pnorm((Y2-Mu2)/Sigma2)

    if (copula == "Gumbel") {

      Phi1 = pmin(0.999,Phi1)
      Phi2 = pmin(0.999,Phi2)

    }

    temp1 = log(c_density(Phi1,Phi2,theta.vec))
    temp2 = ((Y1-Mu1)/Sigma1)^2/2
    temp3 = ((Y2-Mu2)/Sigma2)^2/2
    temp4 = log(2*pi)+log(Sigma1)+log(Sigma2)

    return(temp1-temp2-temp3-temp4)

  }
  normal_LLH = function(para) {

    Mu1 = para[1]
    Mu2 = para[2]

    cov_matrix = array(0,dim = c(2,2,n))
    cov_matrix[1,1,] = Sigma1^2
    cov_matrix[1,2,] = cov_matrix[2,1,] = rho*Sigma1*Sigma2
    cov_matrix[2,2,] = Sigma2^2

    LLH = rep(0,n+1)
    for (i in 1:n) {

      LLH[i] = log(dmvnorm(cbind(Y1[i],Y2[i]),mean = c(Mu1,Mu2),sigma = cov_matrix[,,i]))

    }
    return(LLH)

  }

  ini = c(sum(Y1/Sigma1^2)/sum(1/Sigma1^2),sum(Y2/Sigma2^2)/sum(1/Sigma2^2))

  if (copula != "normal") {

    res = optim(ini,function(para) sum(general_LLH(para)),
                control = list(fnscale = -1),hessian = TRUE)
    LLHV = general_LLH(res$par)
    LLHV[n+1] = sum(LLHV)

  } else {

    res = optim(ini,function(para) sum(normal_LLH(para)),
                control = list(fnscale = -1),hessian = TRUE)
    LLHV = normal_LLH(res$par)
    LLHV[n+1] = sum(LLHV)

    theta.vec = rho
    rho.BC = rho

  }

  V = solve(-res$hessian)

  CommonMean1_res = c(estimate = res$par[1],
                      SE = sqrt(V[1,1]),
                      Lower = res$par[1]-qnorm(1-0.05/2)*sqrt(V[1,1]),
                      Upper = res$par[1]+qnorm(1-0.05/2)*sqrt(V[1,1]))

  CommonMean2_res = c(estimate = res$par[2],
                      SE = sqrt(V[2,2]),
                      Lower = res$par[2]-qnorm(1-0.05/2)*sqrt(V[2,2]),
                      Upper = res$par[2]+qnorm(1-0.05/2)*sqrt(V[2,2]))

  names(LLHV) = c(1:n,"total")

  return(list("Outcome 1" = Y1,
              "Outcome 2" = Y2,
              "Correlation" = rho,
              "Sample size" = n,
              "Copula" = copula,
              "Copula parameter" = theta.vec,
              "Corrected correlation" = rho.BC,
              "CommonMean 1" = CommonMean1_res,
              "CommonMean 2" = CommonMean2_res,V = V,
              "Log-likelihood values" = LLHV))

}
