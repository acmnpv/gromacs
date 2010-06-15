#include "slater_low.h"

cl_R DSlater_1S_1S(cl_R r,cl_R xi,cl_R xj)
{
  cl_R S,rxi,rxj;

  rxi = rxj = S = ZERO;
  rxi = r*xi;
  rxj = r*xj;
  if (xi == xj) {
        if (r == 0LL) {  S = 0LL

    ; } else {  S = -(-33LL*xi + 48LL*exp(2LL*r*xi)*xi - 36LL*r*Power(xi,2LL) - 

          12LL*Power(r,2LL)*Power(xi,3LL))/(24LL*exp(2LL*r*xi)*r) + 

      (-24LL + 24LL*exp(2LL*r*xi) - 33LL*r*xi - 18LL*Power(r,2LL)*Power(xi,2LL) - 

         4LL*Power(r,3LL)*Power(xi,3LL))/(24LL*exp(2LL*r*xi)*Power(r,2LL)) + 

      (xi*(-24LL + 24LL*exp(2LL*r*xi) - 33LL*r*xi - 18LL*Power(r,2LL)*Power(xi,2LL) - 

           4LL*Power(r,3LL)*Power(xi,3LL)))/(12LL*exp(2LL*r*xi)*r)

    ; }
 
  }
  else {
      if (r == 0LL) {  S = 0LL

    ; } else {  S = (exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),3LL) + 

         exp(2LL*r*xj)*Power(xj,4LL)*

          (-3LL*Power(xi,2LL) - r*Power(xi,3LL) + Power(xj,2LL) + r*xi*Power(xj,2LL)) - 

         exp(2LL*r*xi)*Power(xi,4LL)*

          (Power(xi,2LL)*(1LL + r*xj) - Power(xj,2LL)*(3LL + r*xj)))/

       (exp(2LL*r*(xi + xj))*Power(r,2LL)*Power(xi - xj,3LL)*Power(xi + xj,3LL)) + 

      (2LL*(exp(2LL*r*(xi + xj))*Power(Power(xi,2LL) - Power(xj,2LL),3LL) + 

           exp(2LL*r*xj)*Power(xj,4LL)*

            (-3LL*Power(xi,2LL) - r*Power(xi,3LL) + Power(xj,2LL) + r*xi*Power(xj,2LL)) - 

           exp(2LL*r*xi)*Power(xi,4LL)*

            (Power(xi,2LL)*(1LL + r*xj) - Power(xj,2LL)*(3LL + r*xj))))/

       (exp(2LL*r*(xi + xj))*r*Power(xi - xj,3LL)*Power(xi + xj,2LL)) - 

      (2LL*exp(2LL*r*(xi + xj))*(xi + xj)*Power(Power(xi,2LL) - Power(xj,2LL),3LL) + 

         exp(2LL*r*xj)*Power(xj,4LL)*(-Power(xi,3LL) + xi*Power(xj,2LL)) + 

         2LL*exp(2LL*r*xj)*Power(xj,5LL)*

          (-3LL*Power(xi,2LL) - r*Power(xi,3LL) + Power(xj,2LL) + r*xi*Power(xj,2LL)) - 

         exp(2LL*r*xi)*Power(xi,4LL)*(Power(xi,2LL)*xj - Power(xj,3LL)) - 

         2LL*exp(2LL*r*xi)*Power(xi,5LL)*

          (Power(xi,2LL)*(1LL + r*xj) - Power(xj,2LL)*(3LL + r*xj)))/

       (exp(2LL*r*(xi + xj))*r*Power(xi - xj,3LL)*Power(xi + xj,3LL))

    ; }
   
  }
  return S;
}

