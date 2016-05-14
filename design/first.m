g01=load('C:\Users\Alice\Desktop\middle result\evolution\evolution1.txt');
y01=g01(:,1);
y01=log10(y01-2.38095658);
x01=g01(:,3);
y01 = y01(1:50:end, :);
x01 = x01(1:50:end, :);

g02=load('C:\Users\Alice\Desktop\middle result\evolution\evolution2.txt');
y02=g02(:,1);
y02=log10(y02-0.012665233);
x02=g02(:,3);
y02 = y02(1:50:end, :);
x02 = x02(1:50:end, :);

g03=load('C:\Users\Alice\Desktop\middle result\evolution\evolution3.txt');
y03=g03(:,1);
y03=log10(y03-2994.471066);
x03=g03(:,3);
y03 = y03(1:50:end, :);
x03 = x03(1:50:end, :);

g04=load('C:\Users\Alice\Desktop\middle result\evolution\evolution4.txt');
y04=g04(:,1);
y04=log10(y04-263.8958434);
x04=g04(:,3);
y04 = y04(1:50:end, :);
x04 = x04(1:50:end, :);

g05=load('C:\Users\Alice\Desktop\middle result\evolution\evolution5.txt');
y05=g05(:,1);
y05=log10(y05+31025.56024);
x05=g05(:,3);
y05 = y05(1:50:end, :);
x05 = x05(1:50:end, :);

% x=0:0.01:2*pi;
% y1=sin(x);
% y2=cos(x);
% plot(x,y1);
% hold on;
% plot(x,y2,'r');
% axis([0,8,-1.5,1.5]);
% legend('sin(x)','cox(x)');
% text(x(fix(end/2)),y1(fix(end/2)),'\leftarrow sin(x)');
% text(x(fix(end/2)),y2(fix(end/2)),'\leftarrow cos(x)')
%plot(x01,y01,x02,y02,x03,y03,x04,y04);
plot(x01,y01,'bo-',x02,y02,'bx-',x03,y03,'bsquare-',x04,y04,'b+-', x05, y05, 'b*-');
legend('welded beam ','spring','speed reducer','three-bar truss', 'himmelblau');
xlabel('FES'); 
ylabel('log(f(x)-f(x*))');
% text(g01(fix(end/2)),g02(fix(end/2)),'\leftarrow welded beam');
% text(g02(fix(end/2)),g02(fix(end/2)),'\leftarrow spring');
% text(g03(fix(end/2)),g03(fix(end/2)),'\leftarrow speed reducer');
% text(g04(fix(end/2)),g04(fix(end/2)),'\leftarrow three-bar truss');
% text(g05(fix(end/2)),g05(fix(end/2)),'\leftarrow himmelblau');

gtext('welded beam');
gtext('spring');
gtext('speed reducer');
gtext('three-bar truss');
gtext('himmelblau');