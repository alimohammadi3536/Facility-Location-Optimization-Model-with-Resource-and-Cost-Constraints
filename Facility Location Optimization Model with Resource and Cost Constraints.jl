using JuMP, GLPK

# ایجاد مدل بهینه‌سازی
model = Model(GLPK.Optimizer)

# تعریف مجموعه‌ها
I = 1:3
J = 1:15
J1 = 1:8
J2 = 1:6
J3 = 1:3
H = 1:8
T = 1:30

# تعریف پارامترهای ثابت
w1 = 0.5
w2 = 0.5
fsJ1 = [2e8, 5e8]  # مقادیر fsJ1
ftJ2 = [3e8, 6e8]  # مقادیر ftJ2
vs = 30000
vt = 70000
fi = 20000
ksi = 20000
costh = [5e5, 7e5]  # مقادیر costh
costJ3 = [1e9, 3e9]  # مقادیر costJ3
costih = [5e4, 7e4]  # مقادیر costih
costhJ3 = [2e3, 3e3]  # مقادیر costhJ3
costit = [3e5, 4e5]  # مقادیر costit
djj1 = [10, 50]  # مقادیر djj1
djh = [10, 50]  # مقادیر djh
dhj3 = [10, 50]  # مقادیر dhj3
dj1j2 = [10, 50]  # مقادیر dj1j2
d1max = 150
d2max = 200
d3max = 150
d4max = 500
cap1j1 = [3e3, 5e3]  # مقادیر cap1j1
cap2j2 = [3e3, 5e3]  # مقادیر cap2j2
cap3h = [3e3, 5e3]  # مقادیر cap3h
cap4j3 = [3e3, 5e3]  # مقادیر cap4j3
POPij = [1500, 500000]  # مقادیر POPij
Zt = 6
Fbar = 0.5
Fhat = 0
GAMMA = 1
M = 10000000000
n1 = 100
n2 = 80
U1 = 1000
U2 = 800

# تعریف ماتریس‌ها و متغیرها
fit = rand(3, 30)  # ماتریس fit با ابعاد 3x30 و مقادیر بین 0 و 1

# تعریف متغیرهای مدل
@variable(model, GM[J1] >= 0, Bin)
@variable(model, IS[J2] >= 0, Bin)
@variable(model, u[J, J1] >= 0, Bin)
@variable(model, gamma[J1, J2] >= 0, Bin)
@variable(model, D[H] >= 0, Bin)
@variable(model, Lambda[J3] >= 0, Bin)
@variable(model, e[J, H] >= 0, Bin)
@variable(model, f[H, J3] >= 0, Bin)
@variable(model, y[T] >= 0, Bin)
@variable(model, c[I] >= 0)
@variable(model, W[I, J, J1] >= 0)
@variable(model, omega[I, J, H] >= 0)
@variable(model, cs[J1] >= 0, Int)
@variable(model, ct[J2] >= 0, Int)
@variable(model, X1[I] >= 0, Int)
@variable(model, psi[J1, J2] >= 0, Int)
@variable(model, X2[I] >= 0, Int)
@variable(model, PSI[H, J3] >= 0, Int)
@variable(model, X3[I] >= 0, Int)
@variable(model, Dev[T] >= 0, Int)
@variable(model, Ro[I] >= 0)
@variable(model, R1 >= 0)
@variable(model, R2 >= 0)
@variable(model, Pi >= 0)

# تعریف تابع هدف
@objective(model, Min, 
    w1 * (
        sum(fsJ1[mod1(j1, length(fsJ1))] * GM[j1] for j1 in J1) +
        sum(ftJ2[mod1(j2, length(ftJ2))] * IS[j2] for j2 in J2) +
        sum(vs * cs[j1] for j1 in J1) +
        sum(vt * ct[j2] for j2 in J2) +
        sum(fit[1, 2] * djj1[mod1(j, length(djj1))] * X1[i] for j in J for i in I) +
        sum(ksi * dj1j2[mod1(j1 + j2 - 1, length(dj1j2))] * gamma[j1, j2] for j1 in J1 for j2 in J2) +
        sum(costh[mod1(h, length(costh))] * D[h] for h in H) +
        sum(costJ3[mod1(j3, length(costJ3))] * Lambda[j3] for j3 in J3) +
        sum(costih[mod1(i, length(costih))] * X2[i] for i in I) +
        sum(costhJ3[mod1(h + j3 - 1, length(costhJ3))] * PSI[h, j3] for h in H for j3 in J3) +
        sum(costit[mod1(i, length(costit))] * X3[i] for i in I) +
        sum(Dev[t] for t in T)
    ) + w2 * (-sum(c[i] for i in I))
)

# تعریف محدودیت‌ها
@constraint(model, n1 <= sum(X1[i] for i in I) <= U1)
@constraint(model, n2 <= sum(X2[i] for i in I) <= U2)
@constraint(model, [t in T], (sum(fit[1, 2] * X1[i] for i in I) + sum(fit[1, 2] * X3[i] for i in I)) - Zt <= M * y[t])
@constraint(model, [t in T], -(sum(fit[1, 2] * X1[i] for i in I)) + sum(fit[1, 2] * X3[i] for i in I) + Zt <= M * (1 - y[t]))
@constraint(model, [t in T], (sum(fit[1, 2] * X3[i] for i in I)) - (Zt - sum(fit[1, 2] * X1[i] for i in I)) <= Dev[t] + M * (1 - y[t]))
@constraint(model, [i in I], sum(fit[1, 2] * y[t] for t in T) == c[i])
@constraint(model, [j in J], sum(u[j, j1] for j1 in J1) <= 1)
@constraint(model, [j in J], M * sum(u[j, j1] for j1 in J1) >= sum(sum(W[i, j, j1] * POPij[mod1(i + j - 1, length(POPij))] for j1 in J1) for i in I))
@constraint(model, [j in J, h in H], e[j, h] <= D[h])
@constraint(model, [j in J, j1 in J1], u[j, j1] * djj1[mod1(j1, length(djj1))] <= d1max)
@constraint(model, [j in J, h in H], e[j, h] * djh[mod1(h, length(djh))] <= d3max)
@constraint(model, [j1 in J1, j2 in J2], gamma[j1, j2] * dj1j2[mod1(j1 + j2 - 1, length(dj1j2))] <= d2max)
@constraint(model, [h in H, j3 in J3], f[h, j3] * dhj3[mod1(h + j3 - 1, length(dhj3))] <= d4max)
@constraint(model, [j1 in J1], sum(gamma[j1, j2] for j2 in J2) <= 1)
@constraint(model, [h in H], sum(f[h, j3] for j3 in J3) <= 1)
@constraint(model, [j1 in J1], sum(u[j, j1] for j in J) <= GM[j1])
@constraint(model, [h in H], sum(e[j, h] for j in J) >= D[h])
@constraint(model, [j1 in J1], sum(gamma[j1, j2] for j2 in J2) >= GM[j1])
@constraint(model, [h in H], sum(f[h, j3] for j3 in J3) >= D[h])
@constraint(model, [j2 in J2], sum(gamma[j1, j2] for j1 in J1) >= IS[j2])
@constraint(model, [j3 in J3], sum(f[h, j3] for h in H) >= Lambda[j3])
@constraint(model, [i in I], sum(sum(W[i, j, j1] * POPij[mod1(i + j - 1, length(POPij))] for j in J) for j1 in J1) == X1[i])
@constraint(model, [i in I], sum(sum(omega[i, j, h] * POPij[mod1(i + j + h - 1, length(POPij))] for j in J) for h in H) == X2[i])
@constraint(model, sum(Fbar * X2[i] for i in I) + GAMMA * Pi + sum(Ro[i] for i in I) == sum(X3[i] for i in I))
@constraint(model, [i in I], GAMMA + Pi >= Fhat * X2[i])
@constraint(model, [j1 in J1], cs[j1] <= cap1j1[mod1(j1, length(cap1j1))] * GM[j1])
@constraint(model, [h in H], sum(X2[i] for i in I) <= cap3h[mod1(h, length(cap3h))] * D[h])
@constraint(model, [j2 in J2], ct[j2] <= cap2j2[mod1(j2, length(cap2j2))] * IS[j2])
@constraint(model, [j3 in J3], sum(X3[i] for i in I) <= cap4j3[mod1(j3, length(cap4j3))] * Lambda[j3])

# حل مدل
optimize!(model)

# چاپ وضعیت حل و مقادیر
println(termination_status(model))
println("Objective Value: ", objective_value(model))
println("R1 = ", value(R1))
println("R2 = ", value(R2))
println("X1_1 = ", value(X1[1]))
println("X1_2 = ", value(X1[2]))
println("X1_3 = ", value(X1[3]))
println("X2_1 = ", value(X2[1]))
println("X2_2 = ", value(X2[2]))
println("X2_3 = ", value(X2[3]))
