use directories::UserDirs;

fn main() {
    let ud = UserDirs::new().unwrap();
    let home = ud.home_dir();
    let header_file = home
        .join(".local/include/h3/h3api.h")
        .into_os_string()
        .into_string()
        .unwrap();
    let lib_path = home
        .join(".local/lib")
        .into_os_string()
        .into_string()
        .unwrap();
    println!("cargo:rustc-link-search={lib_path}");
    println!("cargo:rustc-link-lib=h3");
    let bindings = bindgen::Builder::default()
        // The input header we would like to generate
        // bindings for.
        .header(header_file)
        // Tell cargo to invalidate the built crate whenever any of the
        // included header files changed.
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        // Finish the builder and generate the bindings.
        .generate()
        // Unwrap the Result and panic on failure.
        .expect("Unable to generate bindings");

    bindings
        .write_to_file("src/libh3_sys.rs")
        .expect("Couldn't write bindings!");
}
