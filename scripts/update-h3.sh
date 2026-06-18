#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(git rev-parse --show-toplevel 2>/dev/null || true)"
if [[ -z "${ROOT_DIR}" ]]; then
  echo "Error: not inside a git repository." >&2
  exit 1
fi
cd "${ROOT_DIR}"

if ! command -v updatecli >/dev/null 2>&1; then
  echo "Error: updatecli not found in PATH." >&2
  exit 1
fi

if ! git diff --cached --quiet --; then
  echo "Error: staged changes detected. Please commit or unstage them first." >&2
  exit 1
fi

extract_package_version_from_stdin() {
  awk '
    BEGIN { in_package = 0 }
    /^\[package\]$/ { in_package = 1; next }
    /^\[/ { if (in_package) exit }
    in_package && $0 ~ /^version = "[0-9]+\.[0-9]+\.[0-9]+"/ {
      line = $0
      sub(/^version = "/, "", line)
      sub(/"$/, "", line)
      print line
      exit
    }
  '
}

extract_h3_version_from_stdin() {
  sed -nE 's/^targz_root_dir = "h3-([0-9]+\.[0-9]+\.[0-9]+)"/\1/p' | head -n1
}

set_package_version() {
  local new_version="$1"
  local tmp_file
  tmp_file="$(mktemp)"
  awk -v new_version="${new_version}" '
    BEGIN { in_package = 0; replaced = 0 }
    /^\[package\]$/ { in_package = 1; print; next }
    /^\[/ { if (in_package) in_package = 0 }
    in_package && !replaced && $0 ~ /^version = "[0-9]+\.[0-9]+\.[0-9]+"/ {
      print "version = \"" new_version "\""
      replaced = 1
      next
    }
    { print }
  ' Cargo.toml > "${tmp_file}"
  mv "${tmp_file}" Cargo.toml
}

parse_semver() {
  local version="$1"
  local major minor patch
  IFS='.' read -r major minor patch <<< "${version}"
  printf '%s %s %s\n' "${major}" "${minor}" "${patch}"
}

head_cargo_toml="$(git show HEAD:Cargo.toml 2>/dev/null || true)"
old_h3_version="$(printf '%s\n' "${head_cargo_toml}" | extract_h3_version_from_stdin || true)"

updatecli pipeline apply --config updatecli.d

new_h3_version="$(extract_h3_version_from_stdin < Cargo.toml || true)"
current_package_version="$(extract_package_version_from_stdin < Cargo.toml || true)"
pre_bump_package_version="${current_package_version}"

if [[ -z "${current_package_version}" ]]; then
  echo "Error: could not read [package].version from Cargo.toml." >&2
  exit 1
fi

if [[ -n "${old_h3_version}" && -n "${new_h3_version}" && "${old_h3_version}" != "${new_h3_version}" ]]; then
  read -r old_h3_major old_h3_minor old_h3_patch <<< "$(parse_semver "${old_h3_version}")"
  read -r new_h3_major new_h3_minor new_h3_patch <<< "$(parse_semver "${new_h3_version}")"
  read -r lib_major lib_minor lib_patch <<< "$(parse_semver "${current_package_version}")"

  bump_kind=""

  if (( new_h3_major > old_h3_major )); then
    bump_kind="major"
    ((lib_major += 1))
    lib_minor=0
    lib_patch=0
  elif (( new_h3_major < old_h3_major )); then
    echo "Error: H3 version downgrade detected (${old_h3_version} -> ${new_h3_version})." >&2
    exit 1
  elif (( new_h3_minor > old_h3_minor )); then
    bump_kind="minor"
    ((lib_minor += 1))
    lib_patch=0
  elif (( new_h3_minor < old_h3_minor )); then
    echo "Error: H3 version downgrade detected (${old_h3_version} -> ${new_h3_version})." >&2
    exit 1
  elif (( new_h3_patch > old_h3_patch )); then
    bump_kind="patch"
    ((lib_patch += 1))
  elif (( new_h3_patch < old_h3_patch )); then
    echo "Error: H3 version downgrade detected (${old_h3_version} -> ${new_h3_version})." >&2
    exit 1
  fi

  if [[ -n "${bump_kind}" ]]; then
    new_package_version="${lib_major}.${lib_minor}.${lib_patch}"
    set_package_version "${new_package_version}"
    echo "Bumped libh3 version (${bump_kind}): ${current_package_version} -> ${new_package_version}"
  fi
fi

post_bump_package_version="$(extract_package_version_from_stdin < Cargo.toml || true)"

files=(Cargo.toml Cargo.lock README.md)
changed=()

for file in "${files[@]}"; do
  if ! git diff --quiet -- "${file}"; then
    changed+=("${file}")
  fi
done

if [[ "${#changed[@]}" -eq 0 ]]; then
  echo "No H3 updates detected. Nothing to commit."
  exit 0
fi

git add -- "${changed[@]}"

h3_version="${new_h3_version:-$(extract_h3_version_from_stdin < Cargo.toml || true)}"
default_message="chore(h3): update H3 sources"
if [[ -n "${h3_version}" ]]; then
  default_message="${default_message} to v${h3_version}"
fi
if [[ -n "${pre_bump_package_version}" && -n "${post_bump_package_version}" ]]; then
  if [[ "${pre_bump_package_version}" != "${post_bump_package_version}" ]]; then
    default_message="${default_message}; libh3 v${pre_bump_package_version} -> v${post_bump_package_version}"
  else
    default_message="${default_message}; libh3 v${post_bump_package_version}"
  fi
fi

commit_message="${1:-${default_message}}"
git commit -m "${commit_message}"

echo "Committed ${#changed[@]} file(s): ${changed[*]}"
