# Maintainer: Brandon Graham Cobb <brandongrahamcobb@ic>
pkgname=lucy
pkgver=120824
pkgrel=1
pkgdesc="A Discord bot using discord.py"
arch=('any')
url="https://github.com/brandongrahamcobb/zip"
license=('GPL')
depends=('python' 'python-virtualenv')
source=("lucy_${pkgver}.tar.gz")
sha256sums=('SKIP')

prepare() {
    cd "$srcdir/${pkgname}_${pkgver}"
    # Any prepare steps go here
}

build() {
    cd "$srcdir/${pkgname}_${pkgver}"
    # Any build steps go here
}

package() {
    cd "$srcdir/${pkgname}_${pkgver}"
    install -Dm755 "lucy.sh" "$pkgdir/usr/bin/lucy.sh"
    install -Dm644 "README.md" "$pkgdir/usr/share/$pkgname/README.md"
    install -Dm644 "LICENSE" "$pkgdir/usr/share/$pkgname/LICENSE"
    cp -r bot "$pkgdir/usr/share/$pkgname/"
    cp -r resources "$pkgdir/usr/share/$pkgname/"
}

