##Nina Haas
##Utilizes Adams-Bashforth Method to approximate input function.
##Since this method requires 3 initial steps, Heum Method is used first

##AB3
##Input: function, initial value, end time, step size, printing step
##Output: vector of y-values and its corresponding vector of times
function [y,t] = AB3(fun,y0,T,DT,IOSTEP)
  n = 2;
  M = T/DT;
  S = floor(T/(IOSTEP*DT))+1;
  y = zeros(length(y0),S);
  t = zeros(1,S);
  y(:,1) = y0;
  
##Heum Method for first two steps 
  y1 = y0 + (DT/2)*(fun(y0 + DT*fun(y0,0),1))+fun(y0,0);
  if(mod(1,IOSTEP) == 0);
    y(:,n) = y1;
    t(n) = DT;
    n = n+1;
  end
  y2 = y1 + (DT/2)*(fun(y1 + DT*fun(y1,DT),2)) + fun(y1,1);
  if(mod(2,IOSTEP) == 0)
    y(n) = y2;
    t(n) = 2*DT;
    n = n+1;
  end

##Adams-Bashforth Methods for all other steps
  for i = 3: M
    y3 = y2 + (DT/12)*(23*fun(y2, (i-1)*DT) - 16*fun(y1,(i-2)*DT) + 5*fun(y0,(i-3)*DT));
    if(mod(i,IOSTEP) == 0)
      y(:,n) = y3;
      t(n) = i*DT;
      n = n+1;
    end
    y0 = y1;
    y1 = y2;
    y2 = y3;
  end    
  