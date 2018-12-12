function [Y] = F(y,t);
  Y(1) = -y(1) + y(2)*y(3);
  Y(2) = -y(2) + (y(3) - 2)*y(1);
  Y(3) = 1 - y(1)*y(2);
  Y = Y';
end