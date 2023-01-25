function [tris] = tris_per_face(N)
tri = zeros(N,N);
k=1;
for i = 1:N
    for j = 1:N+1-i
        tri(j,i)= k;
        k = k + 1;
    end
end

k = (N-1)*(N-1);
tris = zeros(k,3);
k = 1;
for i = 1:N-1
    for j = 1:N-i
        tri1 = [tri(j,i),tri(j+1,i),tri(j,i+1)];
        tris(k,:) = tri1;
        k = k + 1;
    end
end
for i = 2:N-1
    for j = 1:N-i
        tri1 = [tri(j,i),tri(j+1,i-1),tri(j+1,i)];
        tris(k,:) = tri1;
        k = k + 1;
    end
end
end