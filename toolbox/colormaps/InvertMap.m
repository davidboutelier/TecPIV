function InvertMap(inputArg1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

pattern = '.mat';
replacement = '';
InvertName = regexprep(inputArg1,pattern,replacement)
outputName = [InvertName,'-i',pattern]

DATA = load(inputArg1);
RGB = DATA.RGB;
RGB = flipud(RGB);

save(outputName, 'RGB')

end

