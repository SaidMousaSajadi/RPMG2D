function [A_perms , index] = perms_k(A,k)
if k <= factorial(length(A))
    r = 1 ;
    while true
        randindex = randperm(length(A)) ;
        if r == 1
            index(r,:) = randindex ;
        end
        if r > 1
            clear Flag
            for i = 1:size(index,1)
                Flag(i) = any(index(i,:) ~= randindex) ;
            end
            if all(Flag)
                index = [index ; randindex] ;
            end
        end
        r = r + 1 ;
        if size(index,1) == k
            break
        end
    end
    for i=1:size(index,1)
        A_perms(i,:) = A(index(i,:)) ;
    end
else
    error(['Your chosen k is greater than the number of permutations. ' num2str(k) ' is greater of ' num2str(factorial(length(A))) '.'])
end
end