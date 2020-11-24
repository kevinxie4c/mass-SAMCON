# IDSAMCON
Improve SAMCON using inverse dynamics

## Dependency
Install [Open Dynamics Engine](https://www.ode.org/).
Install [Dart](https://dartsim.github.io/) following the [instructions](https://dartsim.github.io/install_dart_on_ubuntu.html).

## Build
```
cmake .
make
```

## Usage
An example of the task file can be seen in https://github.com/kevinxie4c/inverse_dynamics. Edit the task file and then
```
mass.out [task_file] [use_mass] [rounds] [guided_num] [online]
```

## See Also
https://github.com/kevinxie4c/motion-viewer

https://github.com/kevinxie4c/inverse_dynamics
