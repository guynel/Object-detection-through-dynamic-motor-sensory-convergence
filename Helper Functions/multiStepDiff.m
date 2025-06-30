function [x] = multiStepDiff(x,n)
   % n-lagged differences; still return output if input has missing values,
   % with NaNs only where needed 
   if n == 0 
       return
   end
   n = n+1;
   x(:,n:end) = x(:,n:end)-x(:,1:(end-n+1));
   x(:,1:(n-1)) = NaN;
   x = x/(n-1); % N-lagged derivative, corrected for step-size
end