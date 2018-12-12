##Nina Haas
##creates Lagrange polynomials to fit the given domain

##Input: domain vector, specific x value, index value
##Output: Langrange poly as a function of x
function [f] = Lagrange_poly(xi,x,i)
  f = ones(1,length(x));
  for n=1:length(xi)
    if(xi(i)-xi(n) == 0)
      continue
    else
    f = f .* (x-xi(n))/(xi(i)-xi(n));
    end
  end
end