#[cfg(feature = "uber_h3_from_scratch")]
mod uber_h3_from_scratch {
    use sha2::Digest;
    use std::path::{Path, PathBuf};
    fn download(url: &str, dest: &str, sha512: &str) {
        println!("{url} -> {dest}");
        let file_content = reqwest::blocking::get(url)
            .expect("Download failed")
            .bytes()
            .expect("Decoding failed");

        let mut hasher = sha2::Sha512::new();
        hasher.update(&file_content);
        let checksum = hasher.finalize();
        let sha512 = hex::decode(sha512).expect("Checksum not valid HEX!");
        if &checksum[..] != &sha512[..] {
            panic!(
                "Checksum for {} invalid {:?} != {:?}",
                dest,
                checksum.to_vec(),
                sha512
            )
        }

        let mut file = std::fs::File::create(dest).expect("archive could not be created");
        let mut content = std::io::Cursor::new(file_content);
        std::io::copy(&mut content, &mut file).expect("File couldn't be written");
    }

    fn untar(targz: &str, dest: &str) {
        let tar_gz = std::fs::File::open(targz).expect("Archive not found!");
        let tar = flate2::read::GzDecoder::new(tar_gz);
        let mut archive = tar::Archive::new(tar);
        archive.unpack(dest).expect("Unpacking failed");
    }

    struct H3Sources {
        targz_url: String,
        targz_sha512sum: String,
        targz_root_dir: String,
    }

    fn from_toml(toml: &str) -> H3Sources {
        let contents: String = std::fs::read_to_string(toml).expect("File missing: {toml}");
        let table = contents.parse::<toml::Table>().unwrap();

        let h3_sources = table["package"].as_table().expect("Read package failed")["metadata"]
            .as_table()
            .expect("Read metadata failed!")["h3-sources"]
            .as_table()
            .expect("Read h3-sources failed!");

        let targz_url = h3_sources["targz_url"]
            .as_str()
            .expect("targz_url missing!")
            .to_string();

        let targz_sha512sum = h3_sources["targz_sha512sum"]
            .as_str()
            .expect("targz_sha512sum missing!")
            .to_string();
        let targz_root_dir = h3_sources["targz_root_dir"]
            .as_str()
            .expect("targz_root_dir missing!")
            .to_string();

        H3Sources {
            targz_url,
            targz_sha512sum,
            targz_root_dir,
        }
    }

    pub fn build() -> Result<PathBuf, Box<dyn std::error::Error>> {
        let out_dir = std::env::var("OUT_DIR")?;
        let out_dir = Path::new(&out_dir);
        let src_dir = out_dir.join("uber");

        let info = from_toml("Cargo.toml");
        let targz_name = Path::new(&info.targz_url).iter().last().unwrap();

        std::fs::create_dir_all(&src_dir)?;
        let dwnld_dest = out_dir.join(targz_name);
        if !dwnld_dest.exists() {
            download(
                &info.targz_url,
                &dwnld_dest.display().to_string(),
                &info.targz_sha512sum,
            );
        }
        untar(dwnld_dest.to_str().unwrap(), src_dir.to_str().unwrap());

        let build_type = if cfg!(debug_assertions) {
            "Debug"
        } else {
            "Release"
        };
        let dst = cmake::Config::new(src_dir.join(info.targz_root_dir))
            .define("BUILD_TESTING", "OFF")
            .define("BUILD_GENERATORS", "OFF")
            .define("BUILD_BENCHMARKS", "OFF")
            .define("BUILD_FILTERS", "OFF")
            .define("ENABLE_LINTING", "OFF")
            .define("ENABLE_DOCS", "OFF")
            .define("ENABLE_COVERAGE", "OFF")
            .profile(build_type)
            .build();
        Ok(dst)
    }
}

#[cfg(not(feature = "uber_h3_from_scratch"))]
fn from_env() -> std::path::PathBuf {
    let h3_inst_prefix = std::env::var("H3_INSTALL_PREFIX")
        .expect("Environment variable H3_INSTALL_PREFIX must be set");
    std::path::Path::new(&h3_inst_prefix).to_owned()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    #[cfg(not(feature = "uber_h3_from_scratch"))]
    let h3_inst_prefix = from_env();

    #[cfg(feature = "uber_h3_from_scratch")]
    let h3_inst_prefix = uber_h3_from_scratch::build()?;

    let header_file = h3_inst_prefix
        .join("include/h3/h3api.h")
        .display()
        .to_string();
    println!(
        "cargo:rustc-link-search=native={}",
        h3_inst_prefix.join("lib").display()
    );
    println!("cargo:rustc-link-lib=static=h3");
    let bindings = bindgen::Builder::default()
        .header(header_file)
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new())) // invalidate the built crate when header file changed
        .generate()
        .expect("Unable to generate bindings");

    bindings
        .write_to_file("src/libh3_sys.rs")
        .expect("Couldn't write bindings!");
    Ok(())
}
