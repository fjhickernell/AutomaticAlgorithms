%This is a script that runs the orignal foolScript and foolNAGscript
%multiple times with different algorithms and prints to text file that in
%Latex format for easy copying into a paper. 
%The original scripts have actually been turned into functions with a 
%structure variable input and output. Everything runs with the same 
%beginfool, snoopNAG, peakyfunction, and NAGpeaks as the stand alone scripts.

filename='fooltable';
algs=[{'d01ah'},{'d01aj'},{'d01ak'},{'quadgk'},{'quad'},{'chebfun'}];
xlswrite('fooltable',algs,'Sheet1','B1')
categories=[{'Alg'},{'Broken'},{'Real'},{'Error'},...
    {'Infinity Ratio'},{'Two Ratio'}]';
tableinfo.morepeaks=false;
tableinfo.lower=0;
tableinfo.upper=1;
tableinfo.c=100;
tableinfo.p=5;
tableinfo.coefficient=20;
tableinfo.degree=2;
tableinfo.sign=1;
tableinfo.RegFunc=@(x) tableinfo.coefficient.*x.^tableinfo.degree;
tableinfo.RegFuncprime=... %and its derivative
    @(x) tableinfo.degree.*tableinfo.coefficient*x.^(tableinfo.degree-1); %slowly varying
tableinfo.RegFuncdubprime=... %and its second derivative
    @(x) tableinfo.degree.*(tableinfo.degree-1).*tableinfo.coefficient*x.^(tableinfo.degree-2);
% xlswrite('fooltable',categories,'Sheet1','A1')

fID=fopen('tableOfools.txt','w');
fprintf(fID,'\\begin{tabular}{r|c|c|c|c|c}');
formatequ=['\r\n\\multicolumn{6}{c}{\\(f(x)=',num2str(tableinfo.coefficient),'x^',num2str(tableinfo.degree),'\\)}\\\\'];
fprintf(fID,formatequ);
fprintf(fID,'\r\n\\midrule');
fprintf(fID,'\r\nMethod & Failed & Real & RelError & InfinityRatio & TwoRatio \\\\ \r\n\\toprule') ;



tableinfo.quadtype='d01ah';
naginfo=tableNAGscript(tableinfo);
data=[naginfo.secondattempt,naginfo.realintegral,naginfo.error,...
    naginfo.inormratio,naginfo.ratio2norm]';
% xlswrite(filename,data,'Sheet1','B2')
formatspec=['\r\n\\midrule \r\n',tableinfo.quadtype,' & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f \\\\'];
fprintf(fID,formatspec,data);

tableinfo.quadtype='d01aj';
naginfo=tableNAGscript(tableinfo);
data=[naginfo.secondattempt,naginfo.realintegral,naginfo.error,...
    naginfo.inormratio,naginfo.ratio2norm]';
% xlswrite(filename,data,'Sheet1','C2')
formatspec=['\r\n\\midrule \r\n',tableinfo.quadtype,' & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f \\\\'];
fprintf(fID,formatspec,data);

tableinfo.quadtype='d01ak';
naginfo=tableNAGscript(tableinfo);
data=[naginfo.secondattempt,naginfo.realintegral,naginfo.error,...
    naginfo.inormratio,naginfo.ratio2norm]';
% xlswrite(filename,data,'Sheet1','D2')
formatspec=['\r\n\\midrule \r\n',tableinfo.quadtype,' & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f \\\\'];
fprintf(fID,formatspec,data);

tableinfo.fname='quadgk';
info=tableScript(tableinfo);
data=[info.secondattempt,info.realintegral,info.error,...
    info.inormratio,info.ratio2norm]';
% xlswrite(filename,data,'Sheet1','E2')
formatspec=['\r\n\\midrule \r\n',tableinfo.fname,' & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f \\\\'];
fprintf(fID,formatspec,data);

tableinfo.fname='quad';
info=tableScript(tableinfo);
data=[info.secondattempt,info.realintegral,info.error,...
    info.inormratio,info.ratio2norm]';
% xlswrite(filename,data,'Sheet1','F2')
formatspec=['\r\n\\midrule \r\n',tableinfo.fname,' & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f \\\\'];
fprintf(fID,formatspec,data);

tableinfo.fname='chebint';
info=tableScript(tableinfo);
data=[info.secondattempt,info.realintegral,info.error,...
    info.inormratio,info.ratio2norm]';
%xlswrite(filename,data,'Sheet1','G2')
formatspec=['\r\n\\midrule \r\n',tableinfo.fname,' & %8.4f & %8.4f & %8.4f & %8.4f & %8.4f \\\\'];
fprintf(fID,formatspec,data);

fprintf(fID,'\r\n\\end{tabular}');
fclose(fID);