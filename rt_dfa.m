function result=rt_dfa(data,p);
%function result=rt_dfa(data,p);
%Real-time* (line) detrended flux analysis of data with window size 6^p
%
%*Not really real-time though, because we feed in the entire input, but we will simulate as if we would be real-time
%
%    RT_Fractools: Real-time fractal signal processing in the time domain
%    Copyright (C) {2013} Andras Hartmann <andras.hartmann@gmail.com>
%
%    Please, if you use this program for any scientific purpose, cite 
%    Hartmann, A., Mukli, P., Nagy, Z., Kocsis, L., Hermán, P., & Eke, A. (2013).
%    Real-time fractal signal processing in the time domain. Physica A:
%    Statistical Mechanics and Its Applications, 392(1), 89–102.
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

%helper functions just for modularity (keep it simple)
%this is to accumulate the result at the end of any window, not neccessary end of the sliding window
function dfa_helper_0()
    mindex(j) = mindex(j) + 1;
    if (j==1)
	sx(j)=sx(j)+x;
	sx2(j)=sx2(j)+x.*x;
	ssx(j,mindex(j)) = sx(j);
	ssx2(j,mindex(j)) = sx2(j);
    end;
    if (mod(mindex(j),2)==0)
	ssx(j+1, mindex(j)./2) = ssx(j,mindex(j)) + ssx(j,mindex(j)-1);
	ssx2(j+1, mindex(j)./2) = ssx2(j,mindex(j)) + ssx2(j,mindex(j)-1);
    end;
    sxi(j)=sxi(j)+x.*nn;

    % online parameter calculation of line detrending
    m = -6 ./ nn .* (-2 .* sxi(j) + ssx(j, mindex(j)) .* nn + ssx(j, mindex(j))) ./ (nn .^ 2 - 1);
    b = ssx(j, mindex(j)) ./ nn - m .* nn ./ 2 - m ./ 2;

    %The magic again, this one is even more nasty
    ssxd(j,mindex(j)) = (m .^ 2 .* nn .^ 3 ./ 3 + m .^ 2 .* nn .^ 2 ./ 2 + m .^ 2 .* nn ./ 6 + m .* b .* nn .^ 2 + m .* b .* nn + nn .* b .^ 2 + ssx2(j, mindex(j)) - 2 .* sxi(j) .* m - 2 .* ssx(j, mindex(j)) .* b );

%                                                                     /  n                                  \
%                       2  3    2  2    2                             |-----                                |
%                      m  n    m  n    m  n        2              2   | \         2                         |
%                      ----- + ----- + ---- + m a n  + m a n + n a  + |  )   (x(i)  - 2 x(i) m i - 2 x(i) a)|
%                        3       2      6                             | /                                   |
%                                                                     |-----                                |
%                                                                     \i = 1                                /
if ssxd(j,mindex(j))<0
	ssxd(j,mindex(j)) = 0;
	%warning('Results might be compromised because of numerical instability');
end

end;

%this one is for accumulating x and x^2 into sx and sx2
function dfa_helper_1()
    xstart(j)=x;
    if (j==1)
	sx(j)=x;
	sx2(j)=x.*x;
    end;
sxi(j)=x.*nn;
end;

%otherwise: this is the normal functionality accumulating within the window
function dfa_helper_otherwise()
    if (j==1)
	sx(j)=sx(j)+x;
	sx2(j)=sx2(j)+x.*x;
    end;
    sxi(j)=sxi(j)+x.*nn;
end;


%this one calculates the result at the end of the window lgsd is the log of SD 
%FIXME: should not be called SD, sry copy-paste...
function dfa_helper_res()
    for k=1:n;
	lgsd(k)=log(sum(ssxd(k,:))./N);
	lgn(k)=log(2^(k+used1-1));
    end;
    pp=polyfit(lgn,lgsd,1);
    result(value,2) = i;
    result(value,1)=pp(1) ./ 2;
    value = value + 1;
end;

%N=data length
N=2^p;

ignored = [7 1 1; 8 2 1; 9 2 2; 10 2 2; 11 2 2; 12 2 3; 13 2 4; 14 3 5; 15 3 5; 16 3 5; 17 3 6];
% ignored: length, smallest bin size, largest bin size;
i=1;
while (p>ignored(i,1) && i < length(ignored))
    i = i+1;
end;
%FIXME: there should be an error message here if no size found ...

smallest=ignored(i,2);
largest=ignored(i,3);
used1=smallest+1;
used2=p-largest;

%n = number of used points
n=used2-used1+1;

%summed values
sx=zeros(n);
sx2=zeros(n);
sxi=zeros(n);
ssx=zeros(n,n);
ssx2=zeros(n,n);
ssxd=zeros(n,n);
mindex=zeros(n);
value = 1;

%in case you would like to measure the time (don't forget to toc too)
%tic

%initializing data
x = 0;

mindex=zeros(n);
for i=1:length(data);
    x=x+data(i);
    for j=1:n
	c = mod(i,2.^(j+used1-1)); %period
	nn = mod(i-1,2.^(j+used1-1))+1; %n
        switch c;
            case 1
		dfa_helper_1;

            case 0
		dfa_helper_0;
		if (i>=N && j == n)
		    dfa_helper_res;
		end;

	   otherwise
		dfa_helper_otherwise;

	end;
    end;
    if (mod(i,N)==0)
	mindex=zeros(n);
    end;

end;

%in case you would like to measure the time
%toc
%etime(clock,t0)
end
