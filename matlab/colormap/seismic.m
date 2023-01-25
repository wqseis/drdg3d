function map = seismic(m)
%   SEISMIC(M) returns an M-by-3 matrix containing a colormap
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(seismic)
%
% Wenqiang Zhang
% McGill University
% Jan 5 2023
% 
% im=imread('rate_color.png');
% c=squeeze(double(im(fix(end/2),:,:)));
% (or)
% c=squeeze(double(:,im(fix(end/2),:)));

%% Check inputs
narginchk(0,1);

if nargin == 1
    validateattributes(m,{'numeric'},{'numel',1});
end

%% Begin Function
if nargin < 1, m = size(get(gcf,'colormap'),1); end
c=[...
0 0 76;
0 0 82;
0 0 87;
0 0 93;
0 0 98;
0 0 104;
0 0 110;
0 0 115;
0 0 121;
0 0 127;
0 0 132;
0 0 138;
0 0 143;
0 0 149;
0 0 155;
0 0 160;
0 0 166;
0 0 172;
0 0 177;
0 0 183;
0 0 188;
0 0 194;
0 0 200;
0 0 205;
0 0 211;
0 0 217;
0 0 222;
0 0 228;
0 0 233;
0 0 239;
0 0 245;
0 0 250;
2 2 255;
10 10 255;
18 18 255;
26 26 255;
34 34 255;
42 42 255;
50 50 255;
58 58 255;
66 66 255;
74 74 255;
82 82 255;
90 90 255;
98 98 255;
106 106 255;
114 114 255;
122 122 255;
130 130 255;
138 138 255;
146 146 255;
154 154 255;
162 162 255;
170 170 255;
178 178 255;
186 186 255;
194 194 255;
202 202 255;
210 210 255;
218 218 255;
226 226 255;
234 234 255;
242 242 255;
250 250 255;
255 250 250;
255 242 242;
255 234 234;
255 226 226;
255 218 218;
255 210 210;
255 202 202;
255 194 194;
255 186 186;
255 178 178;
255 170 170;
255 162 162;
255 154 154;
255 146 146;
255 138 138;
255 130 130;
255 122 122;
255 114 114;
255 106 106;
255 98 98;
255 90 90;
255 82 82;
255 74 74;
255 66 66;
255 58 58;
255 50 50;
255 42 42;
255 34 34;
255 26 26;
255 18 18;
255 10 10;
255 2 2;
251 0 0;
247 0 0;
243 0 0;
239 0 0;
235 0 0;
231 0 0;
227 0 0;
223 0 0;
219 0 0;
215 0 0;
211 0 0;
207 0 0;
203 0 0;
199 0 0;
195 0 0;
191 0 0;
187 0 0;
183 0 0;
179 0 0;
175 0 0;
171 0 0;
167 0 0;
163 0 0;
159 0 0;
155 0 0;
151 0 0;
147 0 0;
143 0 0;
139 0 0;
135 0 0;
131 0 0;
127 0 0;
];
%... Interpolate get requested size for color table
pp=1:(m-1)/(size(c,1)-1):m;
r=interp1(pp,c(:,1),1:m);
g=interp1(pp,c(:,2),1:m);
b=interp1(pp,c(:,3),1:m);
%... Normalize to range [0,1], and divide again by maximum value
% to correct for round-off errors associated with the interpolation.
map=[r' g' b']/255;
map = map/max(map(:));
