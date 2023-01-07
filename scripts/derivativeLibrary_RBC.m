function [farray, fnameArray] = derivativeLibrary_RBC(f, nosSnapshot, orderMax, degree, fname)
% Input Arguments yName filterType N Nfilter libraryDerivativeOrder 
% thresholdAlpha tolerance T skipSnapshot subtractModel SAVE_DIR

% The function does not return f, hence returnf is always false

addpath('../eqsdiscovery/')

sizeF = size(f);

orderArray = [0 0; 1 0; 0 1; 2 0; 0 2; 1 1; 3 0; 0 3; 2 1; 1 2; ...
    4 0; 0 4; 3 1; 2 2; 1 3; 5 0; 0 5; 4 1; 3 2; 2 3; 1 4; ...
    6 0; 0 6; 5 1; 4 2; 3 3; 2 4; 1 5; ...
    7 0; 0 7; 6 1; 5 2; 4 3; 3 4; 2 5; 1 6; ...
    8 0; 0 8; 7 1; 6 2; 5 3; 4 4; 3 5; 2 6; 1 7];

counter = 0;
for count = 1:length(orderArray)
if sum(orderArray(count,:))>orderMax
    continue;
else
    counter = counter + 1;
    for countT = 1:nosSnapshot
        if count == 1
            tempa = ones(sizeF(1),sizeF(2));
            tempb = ' ';
        else
        [tempa, tempb] = derivative_RBC(...
            f(:,:,countT),[orderArray(count,1),orderArray(count,2)],fname);
        end
        farray(:,:,countT,counter) = tempa;
        fnameArray(counter) = string(tempb);
    end
end
end
end
