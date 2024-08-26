# Maintainer: Brandon Graham Cobb <brandongrahamcobb@icloud.com>
pkgname=lucy
pkgver=260824
pkgrel=1
pkgdesc="A Discord bot using discord.py"
arch=('x86_64')
url="https://github.com/brandongrahamcobb/Documents"
license=('GPL')
depends=('python')
source=("CobbBrandonGraham_${pkgname}_${pkgver}.tar.gz")
sha256sums=('936abf6c9be822ae1945dd1a7fff37738a712731d1d99ce10edea40d49cdb139')

prepare() {
    tar -xzf "CobbBrandonGraham_${pkgname}_${pkgver}.tar.gz"
    # Any additional prepare steps can be added here
}

build() {
    :
}

package() {
    cd "src/CobbBrandonGraham_${pkgname}_${pkgver}"

    # Define the python site-packages directory
    local python_version=$(python -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
    local site_packages_dir="$pkgdir/usr/lib/python${python_version}/site-packages"

    # Create the installation directories
    install -dm755 "$site_packages_dir"
    install -dm755 "$pkgdir/usr/share/doc/$pkgname"
    install -dm755 "$pkgdir/usr/share/$pkgname/resources"

    # Install the Python package
    if [ -d "bot" ]; then
        cp -r bot/* "$site_packages_dir"
    else
        echo "Error: bot directory not found."
        return 1
    fi

    # Install the resources folder
    if [ -d "resources" ]; then
        cp -r resources/* "$pkgdir/usr/share/$pkgname/resources"
    else
        echo "Error: resources directory not found."
        return 1
    fi

    # Install the README and LICENSE files
    install -Dm644 README.md "$pkgdir/usr/share/doc/$pkgname/README"
    install -Dm644 LICENSE "$pkgdir/usr/share/doc/$pkgname/LICENSE"

    # Install the requirements.txt file
    install -Dm644 requirements.txt "$pkgdir/usr/share/doc/$pkgname/requirements.txt"
}
