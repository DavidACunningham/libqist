@echo off

for %%F in (*.tf) DO (
    type %%F | find /v "" > %%F.out
    move %%F.out %%F
)
for %%F in (*.tls.pc) DO (
    type %%F | find /v "" > %%F.out
    move %%F.out %%F
)
exit