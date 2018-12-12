##Nina Haas
##interpolates over each given coordinate point to create one
##complete Lagrange polynomial
##simultaneously uses this polynomial to approximate y at x

##Input: vector of domain, vector of corresponding values, specific x value
##Output: approximated value at x
function [y] = Lagrange_interp(xi,yi,x)
  y = zeros(1,length(x));
  n = length(xi);
  for i=1:n
    y = y + yi(i)*Lagrange_poly(xi,x,i);
  end
end
    