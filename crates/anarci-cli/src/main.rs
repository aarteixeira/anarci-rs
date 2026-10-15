fn main() {
    if let Err(e) = anarci_cli::cli(std::env::args().skip(1).collect()) {
        eprintln!("anarci: error: {e}");
        std::process::exit(1);
    }
}
