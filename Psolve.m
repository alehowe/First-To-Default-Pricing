function [P] = Psolve(BCurve,pi,s,yearfracVector,PPrevious,flag)
% Psolve: computes the survival probability by setting the NPV of the floating 
% leg equal to the NPV of the fixed leg of the CDS given the previous
% probabities

% Inputs:
% Bcurve:               discount curve
% pi:                   recovery rate
% s:                    CDS spreads
% yearfracVector:       year fraction vector between T_i-1 and T_i 
%                       with Eu_30_360 convention
% PPrevious:            vector of previous probabilities 
% flag:                 1: neglecting the accrual term
%                       2: considering the accrual term

% Outputs:
% P:                    update survival probability curve
% index:                first index s.t. the survival probabilities are

% Depending on the flag select the case
switch flag
    
    % Neglect accrual
    case 1 
        % If the probability vector is initialized with an element equal to
        % 1 compute it directly
        if isequal(PPrevious,ones(1,1))

            P=(1-pi)/((s*yearfracVector(end))+(1-pi)); 
        
        % Otherwise compute the new survival probability from the previous
        else

            Num=(1-pi)*sum(BCurve(1:end-1).*(PPrevious(1:end-1)-PPrevious(2:end)))...
            -s*sum(yearfracVector(1:end-1).*BCurve(1:end-1).*PPrevious(2:end))...
                +PPrevious(end)*BCurve(end)*(1-pi);

            Denum=BCurve(end)*(s*yearfracVector(end)-pi+1); 

            P=(Num/Denum);  

        end

     % Consider the accrual
     case 2
        % If the probability vector is initialized with an element equal to
        % 1 compute it directly
        if isequal(PPrevious,ones(1,1))

            P=(1-pi-s*0.5*yearfracVector(end))/(1-pi+s*0.5*yearfracVector(end)); 
        
        % Otherwise compute the new survival probability from the previous
        else

            Num=(1-pi)*sum(BCurve(1:end-1).*(PPrevious(1:end-1)-PPrevious(2:end)))...
                -s*sum(yearfracVector(1:end-1).*BCurve(1:end-1).*PPrevious(2:end))...
                    +0.5*s*sum(yearfracVector(1:end-1).*BCurve(1:end-1).*(PPrevious(1:end-1)-PPrevious(2:end)))...
                        +PPrevious(end)*BCurve(end)*((1-pi)-0.5*s*yearfracVector(end));

            Denum=BCurve(end)*(s*yearfracVector(end)-pi+1-0.5*s*yearfracVector(end)); 

            P=(Num/Denum); 

         end

end
