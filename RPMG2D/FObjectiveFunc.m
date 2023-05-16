function Y = FObjectiveFunc(X,AINumber)
[K,~] = size(AINumber) ;
StrFunc = '' ;
for i = 1:K
  StrFunc = [StrFunc num2str(AINumber(i,1)) '.*exp(-((x-' num2str(AINumber(i,2)) ')./' num2str(AINumber(i,3)) ').^2)'] ;
  if i == 1 && K ~= 1
    StrFunc = [StrFunc ' + '] ;
  elseif i == K
  else
  StrFunc = [StrFunc ' + '] ;
  end
end

StrFunc = strrep (StrFunc, "--", "+") ;
StrFunc = strrep (StrFunc, "-+", "-") ;
StrFunc = strrep (StrFunc, "+-", "-") ;
StrFunc = strrep (StrFunc, "++", "+") ;

F = str2func(['@(x) ' StrFunc]);
Y = F(X) ;
end
