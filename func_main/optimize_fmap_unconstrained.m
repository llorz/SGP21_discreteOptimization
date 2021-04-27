function [x, history, all_C21_iter] = optimize_fmap_unconstrained(C21_ini,metric, Ev2, Ev1)

k1 = size(C21_ini, 1); k2 = size(C21_ini, 2);

switch metric
    case 'ortho'
        eval_fmap = @(C) norm(C*C' - eye(k1), 'fro');
    case 'lap'
        eval_fmap = @(C21)norm(C21*diag(Ev2) - diag(Ev1)*C21,'fro');
end

history.x = [];
history.fval = [];
searchdir = [];



func = @(x) eval_fmap(reshape(x,k1,k2));

x0 = C21_ini(:);
options = optimoptions('fminunc','display','off',  'OutputFcn',@outfun);
x = fminunc(func, x0, options);
num = length(history.fval);
all_C21_iter = cell(num,1);
for k = 1:num
    all_C21_iter{k} = reshape(history.x(:,k), k1, k2);
end


    function stop = outfun(x,optimValues,state)
        stop = false;
        
        switch state
            case 'init'
                hold on
            case 'iter'
                % Concatenate current point and objective function
                % value with history. x must be a row vector.
                history.fval = [history.fval; optimValues.fval];
                history.x = [history.x, x];
                % Concatenate current search direction with
                % searchdir.
                searchdir = [searchdir;...
                    optimValues.searchdirection'];
            case 'done'
                hold off
            otherwise
        end
    end
end