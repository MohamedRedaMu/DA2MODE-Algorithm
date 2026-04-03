% cost_cec2017.m
% Cost-function wrapper for the CEC2017 benchmark suite.
% Returns the optimization error relative to the known global optimum.

function [e] = cost_cec2017(x, fNo)

    %% F2 has been deleted, so skip fucntion ; in CEC2017 it istill have place holder for fNo = 2, so increament it  
    %% Note VI
    % we pass fucntion fNo from 1 to 29 , but actually they are from 1 to  30, but F2 are deleted
    if fNo >=2 
        fNo = fNo + 1 ;
    end

    % Increment the number of function evaluations
    global countFE;
    countFE = countFE + 1 ; % inc function evaluations

    % Define the global minimum values for each function number
    globalMins = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000 ];

    
    % Check if the function number is within the valid range
    if fNo < 1 || fNo > length(globalMins)
        e = Inf;
        fprintf('Function number must be between 1 and %d', length(globalMins));
        return
    end
    
    % Get the global minimum for the given function number
    globalMin = globalMins(fNo);
    
    % Call the CEC2020/2021 function 
    % The solution vector x is transposed to suit the CEC evaluation function requirements
    f = cec17_func(x', fNo); % Standard CEC2017 Library

    % Calculate the error as the difference between the evaluated function value and the global minimum
    e = f - globalMin;
    
end
