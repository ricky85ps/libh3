# libh3 - Wrapper for Uber’s H3 Hexagonal Hierarchical Spatial Index in Rust

This crate calls functions provided by [Uber's H3 library](https://github.com/uber/h3)
to expose a safe Rust API for it.

Contributions are welcome, just do a pull request.

## Quickstart

This refers to a Debian like sytem, please adopt to your needs.
You need [Uber's H3 library](https://github.com/uber/h3) to be installed,
otherwise the libh3 crate cannot link.
To get a version which fits do following steps.
`--prefix $HOME/.local` is where the crate searches for, when building

```bash
sudo apt install git cmake make gcc libtool
git clone --depth=1 --branch v4.2.1 https://github.com/uber/h3.git
cmake -S h3 -B h3_build
cmake --build h3_build --parallel 8
cmake --install ./h3_build --prefix $HOME/.local
```

You may then dir into your project and build it using libh3 as dependency

## Documentation

For further documentation have a look at [https://docs.rs/libh3/](https://docs.rs/libh3/)

For the concepts behind the library refer to [h3geo.org](https://h3geo.org) or
[the H3 blog](https://www.uber.com/en-DE/blog/h3/)
