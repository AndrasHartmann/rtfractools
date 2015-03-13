%Test script for real-time functions

%this is the original H
y(1:2^16) = 0.8;
yy(1:2^16) = 0.2;
y = [y, yy, y];

load signal;

%apply analysis
dfa_res = rt_dfa(x,14);
dfa_res2 = rt_dfa(x,10);

bdswv_res = rt_bdswv(x,14);
bdswv_res2 = rt_bdswv(x,10);

%plotting dfa
figure;
plot(dfa_res(:,2),dfa_res(:,1),'ro',dfa_res2(:,2),dfa_res2(:,1),'g-',dfa_res2(:,2),y(dfa_res2(:,2))','b-');
xlabel('after # data');
ylabel('H');
legend('H estimate with 2^1^4 window size', 'H estimate with 2^1^0 window size', 'true H');
title('Real-Time DFA');

%plotting swv
figure;
plot(bdswv_res(:,2),bdswv_res(:,1),'ro',bdswv_res2(:,2),bdswv_res2(:,1),'g-',bdswv_res2(:,2),y(bdswv_res2(:,2))','b-');
xlabel('after # data');
ylabel('H');
legend('H estimate with 2^1^4 window size', 'H estimate with 2^1^0 window size', 'true H');
title('Real-Time SWV');
