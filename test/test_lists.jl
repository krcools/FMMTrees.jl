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

insert_after!(dl, 4, q)
@test collect(dl) == [1,3,2,4,4]

insert_after!(dl, 100, t)
@test collect(dl) == [1,3,2,100,4,4]

s100 = start(dl); for i in 1:3; _, s100 = next(dl, s100); end
s3 = start(dl); for i in 1:1; _, s3 = next(dl, s3); end
move_before!(dl, s100, s3)
@test collect(dl) == [1,100,3,2,4,4]

sdl = sublist(dl, n, q)
@test collect(sdl) == [3,2]
