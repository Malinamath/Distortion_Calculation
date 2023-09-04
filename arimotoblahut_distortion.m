function arimotoblahut_distortion()

    lamda = [5.0, 2.5, 1.25];
    ERRT = 1e-6;

    for type = 0:2
        fprintf('\nType %d\n', type);
        
        for j = 1:3
            [p, r, q, d] = init_data(type);
            l = lamda(j);

            for i = 1:1000
                error = arimotoblahut_distortion(p, r, q, d, l);
                if error <= ERRT
                    fprintf('Iter: %3d Lamda %5.2lf, Rate = %.5lf, Distortion = %.5lf\n', i, l, calc_rate(p, r, q), calc_dist(p, d, q));
                    break;
                end
            end

            clear p r q d;
        end
    end

    for type = 1:2
        fprintf('\nType %d\n', type);
        count = 0;
        for l = 0:0.001:5
            [p, r, q, d] = init_data(type);
            for i = 1:1000
                error = arimotoblahut_distortion(p, r, q, d, l);
                if error <= ERRT
                    break;
                end
            end
            count = count + i;
            fprintf('Iter: %3d Lamda %5.2lf, Rate = %.5lf, Distortion = %.5lf\n', i, l, calc_rate(p, r, q), calc_dist(p, d, q));
        end
        fprintf('Average iter %.1lf\n', count / (5 / 0.001));
    end

end

function [p, r, q, d] = init_data(type)
    if type == 0
        X = 2;
        XI = 2;
        p = [0.25, 0.75];
        d = [0, 1; 1, 0];
    elseif type == 1
        X = 2;
        XI = 3;
        p = [0.5, 0.5];
        d = [0, 1, 0.25; 1, 0, 0.25];
    elseif type == 2
        X = 2;
        XI = 2;
        p = [0.25, 0.75];
        d = [0, 2; 1, 0];
    end
    
    r = ones(1, XI) / XI;
    q = zeros(XI, X);
end

function error = arimotoblahut_distortion(p, r, q, d, lamda)
    X = length(p);
    XI = length(r);
    
    temp = zeros(1, X);
    c = zeros(1, XI);
    
    for j = 1:X
        for i = 1:XI
            temp(j) = temp(j) + r(i) * exp(-lamda * d(j, i));
        end
    end
    
    for i = 1:XI
        for j = 1:X
            q(i, j) = r(i) * exp(-lamda * d(j, i)) / temp(j);
        end
    end
    
    for i = 1:XI
        for j = 1:X
            r(i) = r(i) + p(j) * q(i, j);
        end
    end
    
    for i = 1:XI
        for j = 1:X
            c(i) = c(i) + p(j) * exp(-lamda * d(j, i)) / temp(j);
        end
    end
    
    error = sum(r .* abs(1 - c));
end

function rate = calc_rate(p, r, q)
    X = length(p);
    XI = length(r);
    
    R = 0;
    for i = 1:XI
        for j = 1:X
            R = R + p(j) * q(i, j) * log(q(i, j) / r(i));
        end
    end
    
    rate = R / log(2.0);
end

function distortion = calc_dist(p, d, q)
    X = length(p);
    XI = size(q, 1);
    
    C = 0;
    for i = 1:XI
        for j = 1:X
            C = C + p(j) * q(i, j) * d(j, i);
        end
    end
    
    distortion = C;
end
