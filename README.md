================================================================================
【 ソフト名 】ChebyshevPolynomialBdG
【  作成者  】@dc1394
================================================================================

★これは何？
　Bogoliubov-de Gennes方程式を数値的に解くコードです。
　このコードは、cometscome様の「ChebyshevPolynomialBdG」（
　https://github.com/cometscome/ChebyshevPolynomialBdG ）を、単純にC++に移植した　ものです。以下は、cometscome様による「ChebyshevPolynomialBdG」のコードの説明で
　す。

--------------------------------------------------------------------------------
　This solves the Bogoliubov-de Gennes equations and gap equations in the
　s-wave superconductor with the use of the Chebyshev polynomial method. See, 
　Y. Nagai, Y. Ota and M. Machida [arXiv:1105.4939 or
　DOI:10.1143/JPSJ.81.024710] 
　http://journals.jps.jp/doi/abs/10.1143/JPSJ.81.024710
--------------------------------------------------------------------------------

　ビルドには、以下のライブラリが必要です。
　・Eigen

★更新履歴
　2016/08/25 ver.0.1   公開。

★ライセンス
　このソフトはフリーソフトウェアです（2条項BSDライセンス）。
--------------------------------------------------------------------------------
　ChebyshevPolynomialBdG
　Copyright (C) 2016 @dc1394

　1.ソースコード形式であれバイナリ形式であれ、変更の有無に関わらず、以下の条件を
　満たす限りにおいて、再配布および利用を許可します。

　1-1.ソースコード形式で再配布する場合、上記著作権表示、 本条件書および第2項の責
　任限定規定を必ず含めてください。
　1-2.バイナリ形式で再配布する場合、上記著作権表示、 本条件書および下記責任限定
　規定を、配布物とともに提供される文書 および/または 他の資料に必ず含めてくださ
　い。

　2.本ソフトウェアは無保証です。自己責任で使用してください。

  Copyright (c) 2016, @dc1394
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
  * Redistributions of source code must retain the above copyright notice, 
    this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright notice, 
    this list of conditions and the following disclaimer in the documentation 
    and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
--------------------------------------------------------------------------------

　ChebyshevPolynomialBdGにはEigen projectによるEigenを使用しています。こちらのラ　イセンスはMPL2になります。
  thomasfermiにはGNU ProjectによるGNU Scientific Libraryを使用しています。こち
　らのライセンスは GNU General Public License になります。
