using FMMTrees
using Base.Test

d = [1,2,3,4]
dl = list(d)

s = start(dl)
_, t = next(dl, s) # t points at 2
_, n = next(dl, t) # n points at 3
_, q = next(dl, n) # q points at 4

move_before!(dl, n, t)
@test collect(dl) == [1,3,2,4]

insert_after!(dl, q, 4)
@test collect(dl) == [1,3,2,4,4]

insert_after!(dl, t, 100)
@test collect(dl) == [1,3,2,100,4,4]

s100 = start(dl); for i in 1:3; _, s100 = next(dl, s100); end
s3 = start(dl); for i in 1:1; _, s3 = next(dl, s3); end
move_before!(dl, s100, s3)
@test collect(dl) == [1,100,3,2,4,4]

sdl = sublist(dl, n, q)
@test collect(sdl) == [3,2]

l = list([1,2,3,4])

a = [1,2,3,4]
l = list(a)
s, i = start(l), 1
while !done(l,s)
    @test l[s] == a[i]
    i += 1
    _, s = next(l,s)
end

s = start(l)
l[s] = 20
@test a[1] == 20

push!(l,21)
@test collect(l) == [20,2,3,4,21]
