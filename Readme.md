# libh3 - Wrapper for Uber’s H3 Hexagonal Hierarchical Spatial Index in Rust

This crate calls functions provided by [Uber's H3 library](https://github.com/uber/h3)
to expose a safe Rust API for it.

Contributions are welcome, just do a pull request.

## Quickstart

This refers to a Debian like sytem, please adopt to your needs.
It will build [Uber's H3 library](https://github.com/uber/h3) during build step,
otherwise the libh3 crate cannot link.

```bash
sudo apt install cmake make gcc libtool
```

You may then dir into your project and build it using libh3 as dependency

## Provide Own H3 Library

To provide your own H3 Library, you must set `H3_INSTALL_PREFIX` as environment variable, so
libh3 can find it's header and static library.
You must then call cargo with feature flag `uber_h3_from_scratch` **disabled**, e.g.:

```bash
cargo build --no-default-features
```

### Provided By Distribution

If you installed H3 via a package manager, it's likely installed to `/usr`.
Set environment accordingly:

```bash
export H3_INSTALL_PREFIX=/usr
```

### Self build

To build from sources do following steps:

```bash
sudo apt install git cmake make gcc libtool
git clone --depth=1 --branch v4.2.1 https://github.com/uber/h3.git
export H3_INSTALL_PREFIX=$PWD/h3_inst
cmake -S h3 -B h3_build
cmake --build h3_build --parallel 8
cmake --install ./h3_build --prefix $H3_INSTALL_PREFIX
```

## Documentation

For further documentation have a look at [https://docs.rs/libh3/](https://docs.rs/libh3/)

For the concepts behind the library refer to [h3geo.org](https://h3geo.org) or
[the H3 blog](https://www.uber.com/en-DE/blog/h3/)
