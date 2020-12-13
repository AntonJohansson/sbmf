# Slightly Beyond Mean-Field -- or something

## Building
```
$ make release
$ make debug
```
for release/debug mode, will produce ```build/sbmf.a```. Dependencies are currently handled in ```third_party``` and will be built with the project.

## Reading the source

The code is structured in a unity-build manner, all public interface is available in ```include/sbmf/sbmf.h```, and ```src/sbmf/sbmf.c``` acts as the main source file.
